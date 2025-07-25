/*
 *  Copyright (c) 2011 Matthew Arsenault
 *
 *  This file is part of Milkway@Home.
 *
 *  Milkway@Home is free software: you may copy, redistribute and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation, either version 3 of the License, or (at your
 *  option) any later version.
 *
 *  This file is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nbody_lua_misc.h"
#include "nbody_lua.h"
#include "nbody_lua_types.h"
#include "milkyway_lua.h"
#include "nbody_lua_util.h"
#include "nbody_check_params.h"
#include "../../third_party/Expanse/src/bfe.h"

int nbParseBFEParams(lua_State* luaSt, NBodyCtx* ctx, int idx, const char* errMsg)
{
    const char* path = NULL;
    lua_getfield(luaSt, idx, "coeffs_file");
    if (!lua_isstring(luaSt, -1)) {
        if (errMsg) mw_lua_perror(luaSt, "%s", errMsg);
        lua_pop(luaSt, 1);
        return 0;
    }
    path = lua_tostring(luaSt, -1);
    ctx->coeffs_file_path = strdup(path);
    lua_pop(luaSt, 1);
    lua_getfield(luaSt, idx, "nmax");
    ctx->nmax = luaL_optinteger(luaSt, -1, 0);
    lua_pop(luaSt, 1);
    lua_getfield(luaSt, idx, "lmax");
    ctx->lmax = luaL_optinteger(luaSt, -1, 0);
    lua_pop(luaSt, 1);
    ctx->bfeModel = bfe_create_from_file(path);
    if (!ctx->bfeModel) return 0;
    ctx->S_nlm = ctx->bfeModel->S_coeffs;
    ctx->T_nlm = ctx->bfeModel->T_coeffs;
    return 1;
}


/* Try to read item at idx as one of the accepted potential types into ctx. */
int nbGetPotentialTyped(lua_State* luaSt, NBodyCtx* ctx, int idx, const char* errMsg)
{
    if (lua_isnoneornil(luaSt, idx))
    {
        ctx->potentialType = EXTERNAL_POTENTIAL_NONE;
        return 0;
    }
    else if (lua_isfunction(luaSt, idx))
    {
        /* Set the potential type. We will be reevaluating the script
         * later on if we're using this, so don't bother getting the
         * closure now. */
        ctx->potentialType = EXTERNAL_POTENTIAL_CUSTOM_LUA;
        return 0;
    }
    else if (lua_istable(luaSt, idx))
    {
        lua_getfield(luaSt, idx, "type");
        const char* potential_type = lua_tostring(luaSt, -1);
        lua_pop(luaSt, 1);
        if (potential_type && strcmp(potential_type, "bfe") == 0)
        {
            if (!nbParseBFEParams(luaSt, ctx, idx, errMsg))
                return 1;
            ctx->potentialType = EXTERNAL_POTENTIAL_BFE;
            return 0;
        }
    }

    else /* The default kind of potential */
    {
        Potential* tmp;

        tmp = expectPotential(luaSt, idx);
        if (!tmp)
        {
            if (errMsg)
            {
                mw_lua_perror(luaSt, "%s", errMsg);
            }

            return 1;
        }

        ctx->potentialType = EXTERNAL_POTENTIAL_DEFAULT;
        ctx->pot = *tmp;

        return checkPotentialConstants(&ctx->pot);
    }
}


