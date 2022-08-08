/*
 * HotFixGlobal_StokesAssembling.h
 *
 *  Created on: Mar 2, 2017
 *      Author: bartsch
 */

#ifndef INCLUDE_ASSEMBLEROUTINES_HOTFIXGLOBAL_ASSEMBLENSE_H_
#define INCLUDE_ASSEMBLEROUTINES_HOTFIXGLOBAL_ASSEMBLENSE_H_

/**
 * This enum class is, as is stated, part of a hotfix. It influences the
 * Galerking Assembling routines for Navier--Stokes, check src/FE/NSE2D_FixPo.C.
 * The problem in these routines is, that they do not distinguish between assembling
 * WITH and WITHOUT the convective part, but always assemble the convective part,
 * too. Sometimes we do not wish to do so, but want to omit the convective term.
 * (Stokes, Upwinding, MDML using Upwinding) This Hotfix is used as follows:
 * in NSE2D_ParMooN.C a variable 'assemble_nse' in global scope with external linkage
 * is declared. Class NSE2D.C sets it to WITH_CONVECTION or WITHOUT_CONVECTION
 * depending on its needs, and several functions in NSE2D_FixPo.C react by assembling
 * or skipping the convective part.
 *
 * This introduces global state, which is discouraged in C++. But since the interface
 * for the functions in NSE2D_FixPo is very rigid, and also the entire
 * assembling process will be reworked in the near future, we went for this
 * (temporary!) shortcut solution.
 */
enum class Hotfixglobal_AssembleNSE{ WITH_CONVECTION, WITHOUT_CONVECTION};

extern Hotfixglobal_AssembleNSE assemble_nse;


#endif /* INCLUDE_ASSEMBLEROUTINES_HOTFIXGLOBAL_ASSEMBLENSE_H_ */
