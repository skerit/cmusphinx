/* -*- c-basic-offset: 4 -*- */
/**
 * @file search_factory.c
 * @author David Huggins Daines
 */

#include <sphinxbase/ckd_alloc.h>
#include <sphinxbase/cmd_ln.h>
#include <sphinxbase/feat.h>
#include <sphinxbase/fe.h>
#include <sphinxbase/err.h>
#include <sphinxbase/strfuncs.h>

#include <multisphinx/cmdln_macro.h>
#include <multisphinx/search_factory.h>
#include <multisphinx/fwdtree_search.h>
#include <multisphinx/fwdflat_search.h>
#include <multisphinx/state_align_search.h>
#include <multisphinx/latgen_search.h>

#include "search_internal.h"

struct search_factory_s
{
    int refcnt;

    int argc;
    char **argv;

    featbuf_t *fb;
    cmd_ln_t *config;
    logmath_t *lmath;
    acmod_t *acmod;
    ngram_model_t *lm;
    dict_t *dict;
    dict2pid_t *d2p;
    glist_t searches;

    /**
     * Used to determine if we need to clone the acmod yet.
     */
    int build_count;
};

static const arg_t ms_args_def[] =
{ MULTISPHINX_OPTIONS, CMDLN_EMPTY_OPTION };

/* Feature and front-end parameters that may be in feat.params */
static const arg_t feat_defn[] =
{ waveform_to_cepstral_command_line_macro(),
    cepstral_to_feature_command_line_macro(),  CMDLN_EMPTY_OPTION };

/* I'm not sure what the portable way to do this is. */
static int file_exists(const char *path)
{
    FILE *tmp;

    tmp = fopen(path, "rb");
    if (tmp)
        fclose(tmp);
    return (tmp != NULL);
}

static void add_file(cmd_ln_t *config, const char *arg, const char *hmmdir,
        const char *file)
{
    char *tmp = string_join(hmmdir, "/", file, NULL);

    if (cmd_ln_str_r(config, arg) == NULL && file_exists(tmp))
        cmd_ln_set_str_r(config, arg, tmp);
    ckd_free(tmp);
}

static void init_defaults(cmd_ln_t *config)
{
    char const *hmmdir, *featparams;

#ifdef MODELDIR
    char const *lmfile, *dictfile;
    /* Set default acoustic and language models. */
    hmmdir = cmd_ln_str_r(config, "-hmm");
    lmfile = cmd_ln_str_r(config, "-lm");
    dictfile = cmd_ln_str_r(config, "-dict");
    if (hmmdir == NULL && hmmdir_exists(MODELDIR "/hmm/en_US/hub4wsj_sc_8k"))
    {
        hmmdir = MODELDIR "/hmm/en_US/hub4wsj_sc_8k";
        cmd_ln_set_str_r(config, "-hmm", hmmdir);
    }
    if (lmfile == NULL && !cmd_ln_str_r(config, "-fsg")
            && !cmd_ln_str_r(config, "-jsgf")
            && file_exists(MODELDIR "/lm/en_US/hub4.5000.DMP"))
    {
        lmfile = MODELDIR "/lm/en_US/hub4.5000.DMP";
        cmd_ln_set_str_r(config, "-lm", lmfile);
    }
    if (dictfile == NULL && file_exists(MODELDIR "/lm/en_US/cmu07a.dic"))
    {
        dictfile = MODELDIR "/lm/en_US/cmu07a.dic";
        cmd_ln_set_str_r(config, "-dict", dictfile);
    }

    /* Expand acoustic and language model filenames relative to installation path. */
    if (hmmdir && !path_is_absolute(hmmdir) && !hmmdir_exists(hmmdir))
    {
        char *tmphmm = string_join(MODELDIR "/hmm/", hmmdir, NULL);
        cmd_ln_set_str_r(config, "-hmm", tmphmm);
        ckd_free(tmphmm);
    }
    if (lmfile && !path_is_absolute(lmfile) && !file_exists(lmfile))
    {
        char *tmplm = string_join(MODELDIR "/lm/", lmfile, NULL);
        cmd_ln_set_str_r(config, "-lm", tmplm);
        ckd_free(tmplm);
    }
    if (dictfile && !path_is_absolute(dictfile) && !file_exists(dictfile))
    {
        char *tmpdict = string_join(MODELDIR "/lm/", dictfile, NULL);
        cmd_ln_set_str_r(config, "-dict", tmpdict);
        ckd_free(tmpdict);
    }
#endif

    /* Get acoustic model filenames and add them to the command-line */
    if ((hmmdir = cmd_ln_str_r(config, "-hmm")) != NULL)
    {
        add_file(config, "-mdef", hmmdir, "mdef");
        add_file(config, "-mean", hmmdir, "means");
        add_file(config, "-var", hmmdir, "variances");
        add_file(config, "-tmat", hmmdir, "transition_matrices");
        add_file(config, "-mixw", hmmdir, "mixture_weights");
        add_file(config, "-sendump", hmmdir, "sendump");
        add_file(config, "-fdict", hmmdir, "noisedict");
        add_file(config, "-lda", hmmdir, "feature_transform");
        add_file(config, "-featparams", hmmdir, "feat.params");
        add_file(config, "-senmgau", hmmdir, "senmgau");
    }
    /* Look for feat.params in acoustic model dir. */
    if ((featparams = cmd_ln_str_r(config, "-featparams")))
    {
        if (cmd_ln_parse_file_r(config, feat_defn, featparams, FALSE) != NULL)
        {
            E_INFO("Parsed model-specific feature parameters from %s\n", featparams);
        }
    }
}

static int search_factory_initialize(search_factory_t *dcf)
{
    if ((dcf->lmath = logmath_init(
            (float64) cmd_ln_float32_r(dcf->config, "-logbase"), 0, FALSE))
            == NULL)
        return -1;
    if ((dcf->fb = featbuf_init(dcf->config)) == NULL)
        return -1;
    if ((dcf->acmod = acmod_init(dcf->config, dcf->lmath, dcf->fb)) == NULL)
        return -1;
    if ((dcf->dict = dict_init(dcf->config, dcf->acmod->mdef)) == NULL)
        return -1;
    if ((dcf->d2p = dict2pid_build(dcf->acmod->mdef, dcf->dict)) == NULL)
        return -1;

    /* If search engines become pluggable then we would scan for them
     * here.  For now we just build the list from the ones that are
     * known inside multisphinx. */
    dcf->searches = glist_add_ptr(dcf->searches,
            (void *) fwdtree_search_query());
    dcf->searches = glist_add_ptr(dcf->searches,
            (void *) fwdflat_search_query());
    dcf->searches
            = glist_add_ptr(dcf->searches, (void *) latgen_search_query());
    dcf->searches = glist_add_ptr(dcf->searches,
            (void *) state_align_search_query());

    return 0;
}

search_factory_t *
search_factory_init_cmdln(cmd_ln_t *config)
{
    search_factory_t *dcf;

    dcf = ckd_calloc(1, sizeof(*dcf));
    dcf->refcnt = 1;

    dcf->config = cmd_ln_retain(config);
    init_defaults(dcf->config);
    if (search_factory_initialize(dcf) < 0)
        goto error_out;

    return dcf;

    error_out: search_factory_free(dcf);
    return NULL;
}

search_factory_t *
search_factory_init_argv(int argc, char const *argv[])
{
    search_factory_t *dcf;
    int i;

    dcf = ckd_calloc(1, sizeof(*dcf));
    dcf->refcnt = 1;

    dcf->argc = argc;
    dcf->argv = ckd_calloc(argc, sizeof(*dcf->argv));
    for (i = 0; i < argc; ++i)
        dcf->argv[i] = ckd_salloc(argv[i]);

    dcf->config
            = cmd_ln_parse_r(NULL, ms_args_def, dcf->argc, dcf->argv, FALSE);
    if (dcf->config == NULL)
        goto error_out;
    init_defaults(dcf->config);
    if (search_factory_initialize(dcf) < 0)
        goto error_out;

    return dcf;

    error_out: search_factory_free(dcf);
    return NULL;
}

search_factory_t *
search_factory_init(char const *key, char const *val, ...)
{
    search_factory_t *sf;
    char const *k, *v, **argv;
    va_list args;
    int argc;

    va_start(args, val);
    argc = 0;
    while ((k = va_arg(args, char const *)) != NULL)
    {
        v = va_arg(args, char const *);
        if (v == NULL)
            E_ABORT("Odd number of arguments passed to search_factory_init(), giving up");
            ++argc;
        }
        va_end(args);

        argv = ckd_calloc(argc * 2 + 2, sizeof(*argv));
        argv[0] = key;
        argv[1] = val;
        va_start(args, val);
        if (key != NULL)
        argc = 2;
        else
        argc = 0;
        while ((k = va_arg(args, char const *))!= NULL)
        {
            v = va_arg(args, char const *);
            argv[argc++] = k;
            argv[argc++] = v;
        }
        va_end(args);

        sf = search_factory_init_argv(argc, argv);
        ckd_free(argv);
        return sf;
    }

search_factory_t *
search_factory_retain(search_factory_t *dcf)
{
    ++dcf->refcnt;
    return dcf;
}

int search_factory_free(search_factory_t *dcf)
{
    int i;
    if (dcf == NULL)
        return 0;
    if (--dcf->refcnt > 0)
        return dcf->refcnt;
    for (i = 0; i < dcf->argc; ++i)
        ckd_free(dcf->argv[i]);
    ckd_free(dcf->argv);
    acmod_free(dcf->acmod);
    cmd_ln_free_r(dcf->config);
    dict2pid_free(dcf->d2p);
    dict_free(dcf->dict);
    featbuf_free(dcf->fb);
    ngram_model_free(dcf->lm);
    logmath_free(dcf->lmath);
    glist_free(dcf->searches);
    ckd_free(dcf);
    return 0;
}

static searchfuncs_t *
search_factory_find(search_factory_t *dcf, char const *name)
{
    gnode_t *gn;

    for (gn = dcf->searches; gn; gn = gnode_next(gn))
    {
        searchfuncs_t *sf = (searchfuncs_t *) gnode_ptr(gn);
        if (0 == strcmp(sf->name, name))
            return sf;
    }
    return NULL;
}

search_t *
search_factory_create_argv(search_factory_t *dcf, search_t *other,
        char const *name, int argc,
        char const *argv[])
{
    /* Find the matching search module. */
    searchfuncs_t *sf;
    cmd_ln_t *config;
    acmod_t *acmod;

    sf = search_factory_find(dcf, name);
    if (sf == NULL)
    {
        E_ERROR("No search module %s found\n", name);
        return NULL;
    }

    /* Copy and override the configuration with extra parameters. */
    /* FIXME: mmory from arg comes from ... ?? */
    if (argc > 0)
        config = cmd_ln_parse_r(cmd_ln_copy(dcf->config), ms_args_def, argc, (char **)argv, FALSE);
    else
        config = cmd_ln_copy(dcf->config);

    /* FIXME: clone dictionary and language model here... */
    if (dcf->build_count > 0)
    acmod = acmod_copy(dcf->acmod);
    else
    acmod = dcf->acmod;
    ++dcf->build_count;

    return (*sf->init)(other, config, acmod, dcf->d2p);
}

search_t *
search_factory_create(search_factory_t *dcf, search_t *other, char const *name, ...)
{
    char const *k, *v, **argv;
    search_t *search;
    va_list args;
    int argc;

    if (name == NULL)
        return search_factory_create_argv(dcf, other, NULL, 0, NULL);
    va_start(args, name);
    argc = 0;
    while ((k = va_arg(args, char *)) != NULL)
    {
        v = va_arg(args, char *);
        if (v == NULL)
            E_ABORT("Odd number of arguments passed to search_factory_create(), giving up");
            ++argc;
        }
        va_end(args);

        argv = ckd_calloc(argc * 2, sizeof(*argv));
        va_start(args, name);
        argc = 0;
        while ((k = va_arg(args, char *))!= NULL)
        {
            v = va_arg(args, char *);
            argv[argc++] = k;
            argv[argc++] = v;
        }
        va_end(args);

        search = search_factory_create_argv(dcf, other, name, argc, argv);
        ckd_free(argv);
        return search;
    }

featbuf_t *
search_factory_featbuf(search_factory_t *dcf)
{
    return dcf->fb;
}

acmod_t *
search_factory_acmod(search_factory_t *dcf)
{
    return dcf->acmod;
}

ngram_model_t *
search_factory_lm(search_factory_t *dcf)
{
    return dcf->lm;
}

dict2pid_t *
search_factory_d2p(search_factory_t *dcf)
{
    return dcf->d2p;
}
