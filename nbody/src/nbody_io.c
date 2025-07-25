/*
 * Copyright (c) 1993, 2001 Joshua E. Barnes, Honolulu, HI.
 * Copyright (c) 2010-2011 Rensselaer Polytechnic Institute.
 * Copyright (c) 2010-2012 Matthew Arsenault
 * Copyright (c) 2016-2018 Siddhartha Shelton
 * This file is part of Milkway@Home.
 *
 * Milkyway@Home is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Milkyway@Home is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Milkyway@Home.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_config.h"

#include "nbody_util.h"
#include "nbody_io.h"
#include "milkyway_util.h"
#include "nbody_coordinates.h"
#include "nbody_mass.h"
#include <string.h>
#include "nbody_defaults.h"
#include "nbody_types.h"
#include "nbody_likelihood.h"

static void nbPrintSimInfoHeader(FILE* f, const NBodyCtx* ctx, const NBodyState* st)
{
    mwvector cmPos;
    mwvector cmVel;
    cmVel = nbCenterOfMom(st);
    if (st->tree.root)
    {
        cmPos = Pos(st->tree.root);
    }
    else
    {
        cmPos = nbCenterOfMass(st);
    }

    fprintf(f,
            "simple_output = %d\n"
            "hasMilkyway  = %d\n"
            "centerOfMass = %f, %f, %f,   centerOfMomentum = %f, %f, %f,\n",
            ctx->SimpleOutput,
            (ctx->potentialType == EXTERNAL_POTENTIAL_DEFAULT),
            X(cmPos), Y(cmPos), Z(cmPos),
            X(cmVel), Y(cmVel), Z(cmVel)
        );
}

static void nbPrintBodyOutputHeader(FILE* f, const NBodyCtx* ctx, mwbool LBavailable)
{
    if (!ctx->SimpleOutput)
    {
        fprintf(f, "# ignore \t id %22s %22s %22s %22s %22s %22s %22s %22s %22s %22s %22s %22s %22s",
                "x", 
                "y",  
                "z",  
                "l",
                "b",
                "r",
                "v_x",
                "v_y",
                "v_z",
                "mass", 
                "v_los",
                "pmra",
                "pmdec"
            );
        if (LBavailable)
        {
            fprintf(f," %22s %22s\n", "Lambda", "Beta");
        }
        else
        {
            fprintf(f,"\n");
        }
    }
    else
    {
        fprintf(f, "# ignore \t id %22s %22s %22s %22s %22s %22s %22s\n",
                "x",
                "y",
                "z",
                "v_x",
                "v_y",
                "v_z",
                "mass"
            );
    }
}

/* output: Print bodies */
int nbOutputBodies(FILE* f, const NBodyCtx* ctx, const NBodyState* st, const NBodyFlags* nbf)
{
    Body* p;
    mwvector lbr;
    mwvector lambdaBetaR;
    
    real vLOS;
    real lambda_val;
    real beta_val;
    real mu_dec;
    real mu_ra;

    /*Get Histogram Parameters for Lambda-Beta Calculation*/
    NBodyLikelihoodMethod method;
    HistogramParams hp;
    NBHistTrig histTrig;
    mwbool LambdaBetaAvailable = FALSE;
    if (nbGetLikelihoodInfo(nbf, &hp, &method) || method == NBODY_INVALID_METHOD)
    {
        mw_printf("Failed to get Histogram Parameters. Not including Lambda-Beta in Output file.\n");
    }
    else
    {
        LambdaBetaAvailable = TRUE;
        nbGetHistTrig(&histTrig, &hp);
    }

    mwbool isLight = FALSE;
    Body* outputTab = st->bestLikelihoodBodyTab;
    if(!ctx->useBestLike || st->bestLikelihood == DEFAULT_WORST_CASE){
        outputTab = st->bodytab;
    }
    const Body* endp = outputTab + st->nbody;

    nbPrintSimInfoHeader(f, ctx, st);
    nbPrintBodyOutputHeader(f, ctx, LambdaBetaAvailable);

    for (p = outputTab; p < endp; p++)
    {
        fprintf(f, "%8d, %8d,", ignoreBody(p), idBody(p));  /* Print if model it belongs to is ignored */
        
        if (!ctx->SimpleOutput)
        {
            lbr = cartesianToLbr(Pos(p), ctx->sunGCDist);
            vLOS = calc_vLOS(Vel(p), Pos(p), ctx->sunGCDist);
            mu_dec = nbVXVYVZtomuDec(Pos(p),Vel(p),ctx->sunVelx,ctx->sunVely,ctx->sunVelz,ctx->sunGCDist,ctx->NGPdec, ctx->NGPra, ctx->lNCP);
            mu_ra = nbVXVYVZtomuRA(Pos(p),Vel(p),ctx->sunVelx,ctx->sunVely,ctx->sunVelz,ctx->sunGCDist,ctx->NGPdec, ctx->NGPra, ctx->lNCP);
            fprintf(f,
                    " %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f",
                    X(Pos(p)), Y(Pos(p)), Z(Pos(p)),
                    L(lbr), B(lbr), R(lbr),
                    X(Vel(p)), Y(Vel(p)), Z(Vel(p)), Mass(p), vLOS, mu_ra, mu_dec);   
            
            if (LambdaBetaAvailable)
            {
                lambdaBetaR = nbXYZToLambdaBeta(&histTrig, Pos(p), ctx->sunGCDist);
                lambda_val = L(lambdaBetaR);
                beta_val = B(lambdaBetaR);
                fprintf(f,", %22.15f, %22.15f\n", lambda_val, beta_val);
            }
            else
            {
                fprintf(f,"\n");
            }
        }
        else
        {
            fprintf(f,
                    " %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f, %22.15f\n",
                    X(Pos(p)), Y(Pos(p)), Z(Pos(p)),
                    X(Vel(p)), Y(Vel(p)), Z(Vel(p)), Mass(p));
        }
    }

    if (fflush(f))
    {
        mwPerror("Body output flush");
        return TRUE;
    }

    return FALSE;
}

int nbWriteBodies(const NBodyCtx* ctx, const NBodyState* st, const NBodyFlags* nbf)
{
    FILE* f;
    int rc = 0;

    if (!nbf->outFileName)
    {
        return 1;
    }

    f = mwOpenResolved(nbf->outFileName, nbf->outputBinary ? "wb" : "w");
    if (!f)
    {
        mw_printf("Failed to open output file '%s'\n", nbf->outFileName);
        return 1;
    }

    if (nbf->outputBinary)
    {
        mw_printf("Binary output unimplemented\n");
        return 1;
    }
    else
    {
        mw_boinc_print(f, "<bodies>\n");
        rc = nbOutputBodies(f, ctx, st, nbf);
        mw_boinc_print(f, "</bodies>\n");
    }

    if (fclose(f) < 0)
    {
        mwPerror("Error closing output file '%s'", nbf->outFileName);
    }

    return rc;
}

