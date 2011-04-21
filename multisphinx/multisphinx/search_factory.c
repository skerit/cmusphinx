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

#include <multisphinx/cmdln_macro.h>
#include <multisphinx/search_factory.h>

struct search_factory_s
{
    int refcnt;

    featbuf_t *fb;
    cmd_ln_t *config;
    logmath_t *lmath;
    acmod_t *acmod;
    ngram_model_t *lm;
    dict_t *dict;
    dict2pid_t *d2p;
    glist_t searches;
};

static const arg_t ps_args_def[] =
{ POCKETSPHINX_OPTIONS, CMDLN_EMPTY_OPTION };

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

static int hmmdir_exists(const char *path)
{
    FILE *tmp;
    char *mdef = string_join(path, "/mdef", NULL);

    tmp = fopen(mdef, "rb");
    if (tmp)
        fclose(tmp);
    ckd_free(mdef);
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
    char const *hmmdir, *lmfile, *dictfile, *featparams;

#ifdef MODELDIR
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

static int search_factory_initialize(search_factory_t *dcf, cmd_ln_t *config)
{
    dcf->config = cmd_ln_retain(config);
    if ((dcf->lmath = logmath_init(
            (float64) cmd_ln_float32_r(config, "-logbase"), 0, FALSE)) == NULL)
        return -1;
    if ((dcf->fb = featbuf_init(config)) == NULL)
        return -1;
    if ((dcf->acmod = acmod_init(config, dcf->lmath, dcf->fb)) == NULL)
        return -1;
    if ((dcf->dict = dict_init(config, dcf->acmod->mdef)) == NULL)
        return -1;
    if ((dcf->d2p = dict2pid_build(dcf->acmod->mdef, dcf->dict)) == NULL)
        return -1;

    /* If search engines become pluggable then we would scan for them
     * here.  For now we just build the list from the ones that are
     * known inside multisphinx. */
    dcf->searches = glist_add_ptr(dcf->searches, fwdtree_search_query());
    dcf->searches = glist_add_ptr(dcf->searches, fwdflat_search_query());
    dcf->searches = glist_add_ptr(dcf->searches, latgen_search_query());

    return 0;
}

config_t *
search_factory_config(char const *key, char const *val, ...)
{
    char const *k, *v, **argv;
    config_t *config;
    va_list args;
    int argc;

    va_start(args, val);
    argc = 0;
    while ((k = va_arg(args, char const *))!= NULL) {
        v = va_arg(args, char const *);
        if (v == NULL) {
            E_ERROR("Odd number of arguments passed to search_factory_config(), giving up");
            return NULL;
        }
        ++argc;
    }
    va_end(args);

    argv = ckd_calloc(argc, sizeof(*argv));
    va_start(args, val);
    argc = 0;
    while ((k = va_arg(args, char const *))!= NULL) {
        argv[argc++] = k;
        v = va_arg(args, char const *);
        argv[argc++] = v;
    }
    va_end(args);

    config = cmd_ln_parse_r(NULL, ps_args_def, argc, argv, FALSE);
    ckd_free(argv);
    return config;
}

search_factory_t *
search_factory_init(config_t *config)
{
    search_factory_t *dcf;

    dcf = ckd_calloc(1, sizeof(*dcf));
    dcf->refcnt = 1;

    if (search_factory_initialize(dcf, config) < 0) {
        search_factory_free(dcf);
        return NULL;
    }

    return dcf;
}

search_factory_t *
search_factory_retain(search_factory_t *dcf)
{
    ++dcf->refcnt;
    return dcf;
}

int search_factory_free(search_factory_t *dcf)
{
    if (dcf == NULL)
        return 0;
    if (--dcf->refcnt > 0)
        return dcf->refcnt;
    ckd_free(dcf);
    return 0;
}

search_t *
search_factory_create(search_factory_t *dcf, char const *name, ...)
{
    search_t *search;

    /* Find the matching search module. */

    /* Copy and override the configuration with extra parameters. */
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
