/* -*- c-basic-offset: 4; indent-tabs-mode: nil -*- */
/* ====================================================================
 * Copyright (c) 2008 Carnegie Mellon University.  All rights 
 * reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer. 
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * This work was supported in part by funding from the Defense Advanced 
 * Research Projects Agency and the National Science Foundation of the 
 * United States of America, and the CMU Sphinx Speech Consortium.
 *
 * THIS SOFTWARE IS PROVIDED BY CARNEGIE MELLON UNIVERSITY ``AS IS'' AND 
 * ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL CARNEGIE MELLON UNIVERSITY
 * NOR ITS EMPLOYEES BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY 
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * ====================================================================
 *
 */
/**
 * @file sphinx_am_fst.cc
 *
 * Build finite-state transducers from Sphinx acoustic model definition files.
 */

#include <fst/fstlib.h>

namespace sb {
#include <cmd_ln.h>
#include <hash_table.h>
#include <err.h>
#include <ckd_alloc.h>
#include "mdef.h"
};

#include "fstprinter.h"

using namespace fst;
using namespace sb;
using sb::int32; /* Override conflicting typedefs */
typedef VectorFst<LogArc> LogVectorFst;
typedef StdArc::StateId StateId;
typedef StdArc::Label Label;

static const arg_t args[] = {
    { "-mdef",
      ARG_STRING,
      NULL,
      "Model definition input file" },
    { "-binfst",
      ARG_STRING,
      NULL,
      "Binary FST output file" },
    { "-txtfst",
      ARG_STRING,
      NULL,
      "Text FST output file" },
    { "-isym",
      ARG_STRING,
      NULL,
      "Input symbol table output file" },
    { "-osym",
      ARG_STRING,
      NULL,
      "Output symbol table output file" },
    { NULL, 0, NULL, NULL }
};

static int
mdef_fst_phone_str(mdef_t * m, int pid, char *buf)
{
    assert(m);
    assert((pid >= m->n_ciphone) && (pid < m->n_phone));

    buf[0] = '\0';
    sprintf(buf, "%s-%s+%s",
                mdef_ciphone_str(m, m->phone[pid].lc),
                mdef_ciphone_str(m, m->phone[pid].ci),
                mdef_ciphone_str(m, m->phone[pid].rc));
    return 0;
}


static StdVectorFst *
mdef_to_fst(mdef_t *mdef)
{
    StdVectorFst *model = new StdVectorFst;
    SymbolTable *isym = new SymbolTable("triphones");
    SymbolTable *osym = new SymbolTable("phones");
    StateId start = model->AddState();
    StateId end = model->AddState();
    model->SetStart(start);
    model->SetFinal(end, 0);

    int32 offset = mdef_n_ciphone(mdef);
    for (int32 i = offset; i < mdef_n_phone(mdef); ++i) {
	char buf[40];
	mdef_fst_phone_str(mdef, i, buf);
	if (isym->Find(buf) < 0)
            isym->AddSymbol(buf, i-offset+1);
    }    
    
    for (int32 i = 0; i < mdef_n_ciphone(mdef); ++i) {
        osym->AddSymbol(mdef_ciphone_str(mdef, i), i+1);
    }

    model->SetInputSymbols(isym);
    model->SetOutputSymbols(osym);
    return model;
}

int
main(int argc, char *argv[])
{
    mdef_t *mdef = NULL;
    cmd_ln_t *config;

    if ((config = cmd_ln_parse_r(NULL, args, argc, argv, TRUE)) == NULL) {
        return 1;
    }
    if (cmd_ln_str_r(config, "-mdef")) {
        if ((mdef = mdef_init(cmd_ln_str_r(config, "-mdef"), TRUE)) == NULL) {
            E_ERROR("Failed to load model definition file from %s\n",
                    cmd_ln_str_r(config, "-mdef"));
            return 1;
        }
    }

    /* Create FST from dict. */
    StdVectorFst *model = mdef_to_fst(mdef);

    /* Write the model in the requested format. */
    char const *outfile;
    if ((outfile = cmd_ln_str_r(config, "-binfst")) != NULL) {
        model->Write(outfile);
    }
    if ((outfile = cmd_ln_str_r(config, "-txtfst")) != NULL) {
        FstPrinter<StdArc> printer(*model,
                                   model->InputSymbols(),
                                   model->OutputSymbols(),
                                   NULL, false);
        ostream *ostrm = new ofstream(outfile);
        printer.Print(ostrm, outfile);
    }
    if ((outfile = cmd_ln_str_r(config, "-isym")) != NULL) {
        model->InputSymbols()->WriteText(outfile);
    }
    if ((outfile = cmd_ln_str_r(config, "-osym")) != NULL) {
        model->OutputSymbols()->WriteText(outfile);
    }

    return 0;
}
