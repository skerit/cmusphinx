/*
 * search_factory.h
 *
 *  Created on: 2011-03-20
 *      Author: dhuggins
 */

#ifndef __search_FACTORY_H__
#define __search_FACTORY_H__

#include <stdarg.h>

#include <sphinxbase/ngram_model.h>
#include <sphinxbase/cmd_ln.h>

#include <multisphinx/win32_export.h>
#include <multisphinx/search.h>
#include <multisphinx/dict2pid.h>
#include <multisphinx/acmod.h>

typedef struct search_factory_s search_factory_t;

/**
 * Construct a search factory from a list of key-value pairs.
 *
 * This will initialize the various models needed for search initialization.
 */
search_factory_t *search_factory_init(char const *key, char const *val, ...);

/**
 * Construct a search factory from an array of strings.
 */
search_factory_t *search_factory_init_argv(int argc, char const *argv[]);

/**
 * Construct a search factory from a previously parsed command line.
 */
search_factory_t *search_factory_init_cmdln(cmd_ln_t *config);

/**
 * Retain a reference to a search factory.
 */
search_factory_t *search_factory_retain(search_factory_t *dcf);

/**
 * Release a reference to a search factory.
 */
int search_factory_free(search_factory_t *dcf);

/**
 * Construct a new search from the factory, optionally overriding some parameters.
 */
search_t *search_factory_create(search_factory_t *dcf, search_t *other, char const *name, ...);

/**
 * Construct a new search from the factory, optionally overriding some parameters.
 */
search_t *search_factory_create_argv(search_factory_t *dcf, search_t *other, char const *name,
				     int argc, char const *argv[]);

/**
 * Get the feature buffer from the search factory.
 */
featbuf_t *search_factory_featbuf(search_factory_t *dcf);

/**
 * Get the acoustic model from the search factory.
 */
acmod_t *search_factory_acmod(search_factory_t *dcf);

/**
 * Get the language model from the search factory.
 */
ngram_model_t *search_factory_lm(search_factory_t *dcf);

/**
 * Get the dictionary/tied state mapping from the search factory.
 */
dict2pid_t *search_factory_d2p(search_factory_t *dcf);

#endif /* __search_FACTORY_H__ */
