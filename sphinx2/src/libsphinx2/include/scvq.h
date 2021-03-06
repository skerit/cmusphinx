/* ====================================================================
 * Copyright (c) 1999-2004 Carnegie Mellon University.  All rights
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

/*
 * scvq.h -- Interface to semi-continous quantization
 * 
 * **********************************************
 * CMU ARPA Speech Project
 *
 * Copyright (c) 1996 Carnegie Mellon University.
 * ALL RIGHTS RESERVED.
 * **********************************************
 * 
 * HISTORY
 * 
 * $Log$
 * Revision 1.5  2004/12/10  16:48:58  rkm
 * Added continuous density acoustic model handling
 * 
 * 
 * 24-Nov-2004	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University
 * 		Changed return type for SCVQScores and SCVQScores_all to void,
 * 		and removed scvq_sen_psen() and scvq_set_bestpscr() definitions.
 * 		These are now in the senscr module, for integrating with S3
 * 		(continuous) acoustic model handling.
 * 
 * 16-May-96	M K Ravishankar (rkm@cs.cmu.edu) at Carnegie Mellon University
 * 		Created from Fil Alleva's original.
 */

#ifndef __SCVQ_H__
#define __SCVQ_H__

#define CEP_VECLEN	13
#define POW_VECLEN	3
#define CEP_SIZE	CEP_VECLEN
#define POW_SIZE	POW_VECLEN

#define NUM_FEATURES	4

typedef enum {CEP_FEAT=0, DCEP_FEAT=1, POW_FEAT=2, DDCEP_FEAT=3} feat_t;

typedef enum {AGC_NONE=0, AGC_BETA=1, AGC_NOISE=2, AGC_EMAX=3, AGC_MAX=4} scvq_agc_t;
typedef enum {NORM_NONE=0, NORM_UTT=1, NORM_PRIOR=2} scvq_norm_t;
typedef enum {COMPRESS_NONE=0, COMPRESS_UTT=1, COMPRESS_PRIOR=2} scvq_compress_t;

#define MAX_TOPN	6	/* max number of TopN codewords */

#define NONE		-1

void SCVQInit(int32 top, int32 numModels, int32 numDist, double vFloor, int32 use20msdp);
void SCVQNewUtt(void);
void SCVQEndUtt ( void );
void SCVQAgcSet(scvq_agc_t agc_type);
void SCVQAgcInit(int32 agcColdInit, int32 agcLookAhead);
void SCVQSetSenoneCompression (int32 size);
int32 SCVQInitFeat(feat_t feat, char *meanPath, char *varPath, int32 *opdf);
void SCVQSetdcep80msWeight (double arg);
int SCVQComputeFeatures(float **cep,
			float **dcep,
			float **dcep_80ms,
			float **pow,
			float **ddcep,
			float *in);

/* used to fill input buffer */
int32 SCVQNewFrame(float *in);

/* used by search to compute scores */
int32 SCVQGetNextScores(int32 *scores);
void SCVQScores (int32 *scores,
		 float *cep,
		 float *dcep,
		 float *dcep_80ms,
		 float *pcep, float *ddcep);
void SCVQScores_all (int32 *scores,
		     float *cep,
		     float *dcep,
		     float *dcep_80ms,
		     float *pcep,
		     float *ddcep);

void setVarFloor(double aVal);
int32  readMeanCBFile(feat_t feat, float **CB, char *MeanCBFile);
int32  readVarCBFile(feat_t feat, register int32 *det, float **CB, char *VarCBFile);
int32  setPowVar(register int32 *det, float **CB, double pow_var);


#ifndef FALSE
#define FALSE		0
#endif
#ifndef TRUE
#define TRUE		1
#endif

#endif
