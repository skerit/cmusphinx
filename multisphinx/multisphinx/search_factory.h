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

#include <multisphinx/win32_export.h>
#include <multisphinx/search.h>
#include <multisphinx/dict2pid.h>
#include <multisphinx/acmod.h>

typedef struct search_factory_s search_factory_t;
typedef struct cmd_ln_s config_t;

/**
 * Construct a set of configuration parameters.
 */
config_t *search_factory_config(char const *key, char const *val, ...);

/**
 * Construct a search factory from a set of parameters.
 *
 * This will initialize the various models needed for search initialization.
 */
search_factory_t *search_factory_init(config_t *config);

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
search_t *search_factory_create(search_factory_t *dcf, char const *name, ...);

/**
 * Get the feature buffer from the search factory.
 */
featbuf_t *search_factory_fb(search_factory_t *dcf);

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
