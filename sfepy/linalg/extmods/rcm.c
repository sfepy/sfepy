#include "rcm.h"

/*!
  @file
  Reverse Cuthill-McKee algorithm.
*/

/*!
  @fn void rcm_rootls (int32 root, int32 xadj[], int32 adjncy[], int32 mask[],
  int32 *nlvl, int32 xls[], int32 ls[])
  @short Generate a level structure.

  @param root [I] - the node at which the level structure is to be rooted.  
  @param xadj, adjncy [I] - adjacency structure pair for the given graph. 
  @param mask [I] - is used to specify a section subgraph. nodes with
  mask(i)=0 are ignored.                               
  @param nlvl [O] - is the number of levels in the level structure.
  @param xls, ls [O] - array pair for the rooted level structure.           

  rootls() generates the level structure rooted at the input node called
  root. Only those nodes for which mask is nonzero will be considered.

  Adopted from code by Ales Janka.
  
  @par Revision history:
  - 11.10.2000, 0.27.0, c
  - 12.10.2000, 0.27.0
*/
void rcm_rootls (int32 root, int32 xadj[], int32 adjncy[], int32 mask[],
                 int32 *nlvl, int32 xls[], int32 ls[])
{         
  int32 i, j;
  int32 jstop, jstrt;
  int32 lbegin, ccsize, lvlend, lvsize;
  int32 nbr, node;

/*   ------------------                                               
 *   initialization ...                                               
 *   ------------------ */                                              
  mask[root] = 0;                                               
  ls[0] = root;                                                     
  (*nlvl) = 0;                                                         
  lvlend = 0;
  ccsize = 1;

/*  -----------------------------------------------------            
 *  lbegin is the point32er to the beginning of the current            
 *  level, and lvlend point32s to the end of this level.               
 *  ----------------------------------------------------- */           
  do {
    lbegin = lvlend;
    lvlend = ccsize;
    xls[*nlvl] = lbegin;
    (*nlvl)++;
/*  -------------------------------------------------                
 *  generate the next level by finding all the masked                
 *  neighbors of nodes in the current level.                         
 *  ------------------------------------------------- */               
//      printf( "%d %d\n", lbegin, lvlend - 1 );
    for (i = lbegin; i < lvlend; i++) {
      node = ls[i];
      jstrt = xadj[node];
      jstop = xadj[node+1];
//        printf( "  %d %d %d\n", node, jstrt, jstop-1 );
      for (j = jstrt; j < jstop; j++) {
        nbr = adjncy[j];
//  	printf( "     %d %d\n", nbr, mask[nbr] );
        if (mask[nbr] == 0) continue;
	ls[ccsize] = nbr;
        ccsize++;
        mask[nbr] = 0;
      }
    }
           
/*   ------------------------------------------                       
 *   compute the current level width.                                 
 *   if it is nonzero, generate the next level.                       
 *   ------------------------------------------   */                    
    lvsize = ccsize - lvlend;
//      printf( "ss %d\n", lvsize );
  } while (lvsize > 0);
/*   -------------------------------------------------------          
 *   reset mask to one for the nodes in the level structure.          
 *   ------------------------------------------------------- */         
  xls[*nlvl] = lvlend;
//    printf( "dd %d %d\n", *nlvl, xls[*nlvl] );
  for (i = 0; i < ccsize; i++) {
    node = ls[i];
    mask[node] = 1;
  }
  return;
}

/*!
  @fn rcm_fnroot(int32 *root, int32 xadj[], int32 adjncy[], int32 mask[], 
  int32 *nlvl, int32 xls[], int32 ls[])
  @short Find pseudo-peripheral nodes.

  @param xadj, adjncy [I] - adjacency structure pair for the graph.          
  @param mask [I] - specifies a section subgraph. nodes for which              
  mask is zero are ignored by fnroot.                        
  @param root [U] - on input, it (along with mask) defines the                 
  component for which a pseudo-peripheral node is            
  to be found. on output, it is the node obtained.           
  @param nlvl [O] - is the number of levels in the level structure             
  rooted at the node root.                                   
  @param xls, ls [O] - the level structure array pair containing              
  the level structure found.                             

  rcm_fnroot() implements a modified version of the scheme by gibbs, poole, and
  stockmeyer to find pseudo-peripheral nodes. It determines such a node
  for the section subgraph specified by mask and root.

  @par Subroutines:
  rcm_rootls().                       

  Adopted from code by Ales Janka.
  
  @par Revision history:
  - 11.10.2000, 0.27.0, c
*/
void rcm_fnroot (int32 *root, int32 xadj[], int32 adjncy[], int32 mask[],
		 int32 *nlvl, int32 xls[], int32 ls[])
{
  int32 ccsize, j, jstrt, k, kstop, kstrt;
  int32 mindeg, nabor, ndeg, node, nunlvl;

/*  ---------------------------------------------                    
 *  determine the level structure rooted at root.                    
 *  --------------------------------------------- */                   
  rcm_rootls(*root, xadj, adjncy, mask, nlvl, xls, ls);
  ccsize = xls[*nlvl];
  if ((*nlvl == 1) || (*nlvl == ccsize)) return;
/*  ----------------------------------------------------             
 *  pick a node with minimum degree from the last level.             
 *  ---------------------------------------------------- */            
  do {
    jstrt = xls[*nlvl-1];
    mindeg = ccsize;
    (*root) = ls[jstrt];
    if (ccsize != jstrt) {
      for (j = jstrt; j < ccsize; j++) {
        node = ls[j];
        ndeg = 0;
        kstrt = xadj[node];
        kstop = xadj[node+1];
        for (k = kstrt; k < kstop; k++) {
           nabor = adjncy[k];
           if (mask[nabor] > 0)  ndeg++;
        }
        if (ndeg < mindeg) {
          (*root) = node;
           mindeg = ndeg;
        }
      }
    }
/*  ----------------------------------------                         
 *  and generate its rooted level structure.                         
 *  ---------------------------------------- */           
     rcm_rootls (*root, xadj, adjncy, mask, &nunlvl, xls, ls);
     if (nunlvl <= *nlvl)  return;
     *nlvl = nunlvl;
  } while (*nlvl < ccsize);
  return;
}

/*!
  @fn void rcm_rcm (int32 root, int32 xadj[], int32 adjncy[], int32 mask[], int32 perm[],
  int32 ccsize, int32 deg[])
  @short Number a connected component.

  @param root [I] - is the node that defines the connected component and it
  is used as the starting node for the rcm ordering.
  @param xadj, adjncy [I] - adjacency structure pair for the graph.          
  @param deg [I] - is a temporary vector used to hold the degree of the
  nodes in the section graph specified by mask and root.
  @param mask [U] - only those nodes with nonzero input mask values are
  considered by the routine. The nodes numbered by rcm_rcm()
  will have their mask values set to zero.
  @param perm [O] - will contain the rcm ordering.                            
  @param ccsize [O] - is the size of the connected component that has been
  numbered by rcm().                       

  Number a connected component specified by mask and root, using the
  rcm algorithm. The numbering is to be started at the node root.          

  Adopted from code by Ales Janka.
  
  @par Revision history:
  - 11.10.2000, 0.27.0, c
  - 12.10.2000, 0.27.0
*/
void rcm_rcm (int32 root, int32 xadj[], int32 adjncy[], int32 mask[], int32 perm[],
              int32 ccsize, int32 deg[])
{                               
  int32 fnbr, i, j, jstop, jstrt, k, l, lbegin, lnbr, lperm, lvlend;
  int32 nbr, node;

  mask[root] = 0;
  if ( ccsize <= 1 ) return;
  lvlend = 0;
  lnbr = 0;
/*  --------------------------------------------                     
 *  lbegin and lvlend point32 to the beginning and                     
 *  the end of the current level respectively.                       
 *  --------------------------------------------  */
  do {
    lbegin = lvlend;
    lvlend = lnbr + 1;
    for (i = lbegin; i < lvlend; i++) {
/*  ----------------------------------                            
 *  for each node in current level ...                            
 *  ---------------------------------- */
      node = perm[i];
      jstrt = xadj[node];
      jstop = xadj[node+1];
/*  ------------------------------------------------              
 *  find the unnumbered neighbors of node.                        
 *  fnbr and lnbr point32 to the first and last                     
 *  unnumbered neighbors respectively of the current              
 *  node in perm.                                                 
 *  ------------------------------------------------ */ 
      fnbr = lnbr + 1;
      for (j = jstrt; j < jstop; j++) {
	nbr = adjncy[j];
	if ( mask[nbr] != 0 ) {
	  lnbr++;
	  mask[nbr] = 0;
	  perm[lnbr] = nbr;
	}
      }
      if ( fnbr < lnbr ) {
/*  ------------------------------------------                 
 *  sort the neighbors of node in increasing                   
 *  order by degree. linear insertion is used.                 
 *  ------------------------------------------ */
	k = fnbr;
	do {
	  l = k;
	  k++;
	  nbr = perm[k];
	  while ( l >= fnbr ) {
	    lperm = perm[l];
	    if (deg[lperm] <= deg[nbr]) break;
	    perm[l+1] = lperm;
	    l--;
	  }                                        
	  perm[l+1] = nbr;
	} while (k < lnbr);
      }
    }
  } while (lnbr >= lvlend);
/*  ---------------------------------------                         
 *  we now have the cuthill mckee ordering.                         
 *  reverse it below ...                                            
 *  ---------------------------------------  */  

  k = ccsize/2;
  l = ccsize-1;
  for (i = 0; i < k; i++) {
    lperm = perm[l];
    perm[l] = perm[i];
    perm[i] = lperm;
    l--;
  }
  return;
}

#undef __FUNC__
#define __FUNC__ "rcm_genrcm"
/*!
  @fn void rcm_genrcm (int32 neqns, int32 xadj[], int32 adjncy[], int32 perm[])
  @short General reversed Cuthill-McKee permutation algorithm.

  @param neqns [I] - number of equations
  @param xadj, adjncy [I] - array pair containing the adjacency
  structure of the graph of the matrix.
  @param perm [O] - vector that contains the RCM ordering. 
  perm[i] = # means that old #-th node is the new i-th node

  @par Working data:
  \a mask is used to mark variables that have been                  
  numbered during the ordering process. it is               
  initialized to 1, and set to zero as each node            
  is numbered.
  @par                                              
  \a xls is the index vector for a level structure.  the               
  level structure is stored in the currently                
  unused spaces in the permutation vector perm[].             
                                                                          
  General reversed Cuthill-McKee permutation algorithm, general = for a
  graph with more components.
  rcm_genrcm() finds the reverse cuthill-mckee ordering for a general
  graph. For each connected component in the graph, rcm_genrcm obtains
  the ordering by calling the subroutine rcm_rcm().

  \b IMPORTANT \b PROPERTY: If A is block diagonal, then the ordering of the
  blocks is kept unchanged, ie. the numbers of each separate component of
  the graph remain unchanged.

  @par Subroutines:
  rcm_rcm(), rcm_fnroot()

  Adopted from code by Ales Janka.
  
  @par Revision history:
  - 11.10.2000, 0.27.0, c
  - 16.10.2000, 0.27.1
  - 04.06.2001, adopted for rcfem2
*/
void rcm_genrcm(int32 *perm, int32 neqns, int32 *xadj, int32 n_ptr,
		int32 *adjncy, int32 n_indx)
{                                                       
  int32 ccsize, i, nlvl, num, root;
  int32 *xls, *deg, *mask;

  deg = alloc_mem( int32, neqns );
  mask = alloc_mem( int32, neqns );
  xls = alloc_mem( int32, neqns+1 );      /* majorant length*/

  for (i = 0; i < neqns; i++) {
    mask[i] = 1;                    /* no node has been processed yet */
    deg[i] = xadj[i+1] - xadj[i] - 1;      /* degree of the i-th node */
//      printf( "%d %d\n", i, deg[i] );
  }
  num = 0;

  for (i = 0; i < neqns; i++) {
/*  ---------------------------------------                       
 *  for each masked connected component ...                       
 *  ---------------------------------------  */    
    if (mask[i] == 0) continue;                 /* take next (i) */       
    root = i;
/*  -----------------------------------------                  
 *  first find a pseudo-peripheral node root.                  
 *  note that the level structure found by                     
 *  fnroot is stored starting at perm[num].                    
 *  then rcm is called to order the component                  
 *  using root as the starting node.                           
 *  -----------------------------------------  */
    rcm_fnroot (&root, xadj, adjncy, mask, &nlvl, xls, &perm[num]);
    ccsize = xls[nlvl];
//      printf( "%d %d\n", root, ccsize );

    rcm_rcm(root, xadj, adjncy, mask, &perm[num], ccsize, deg);
    num += ccsize;
    if (num > neqns) {
      goto end_label;
    }
  }
 end_label:
  free_mem( deg );
  free_mem( xls );
  free_mem( mask );
  return;
}


#undef __FUNC__
#define __FUNC__ "gr_permuteInPlace"
/*!
  @fn int32 gr_permuteInPlace( GraphCSR *obj, int32 *perm, int32 *permI )
  @short Graph adjacency matrix permutation in place.
  
  @param obj [U] - is the graph to be permuted.
  Some fields are @b overwritten.
  @param perm [I] - is the permutation vector.
  perm[#] = $ means that $ is the old position and # is the new one.
  @param permI [I] - is the inverse permutation vector.

  Permute the adjacency matrix of the graph
  in the place of the original graph structure.

  @par Revision history:
  - 17.10.2000, 0.27.2, c
  - 27.05.2001
*/
int32 gr_permuteInPlace(int32 *row, int32 n_row,
			int32 *col, int32 n_col,
			int32 *perm, int32 n_perm,
			int32 *permI, int32 n_permI)
{
  int32 ret = 0;
  int32 ir, ic, iold, iord;
  int32 tmpcol, neword, auxint;
  int32 nNod = n_perm, nEdge = n_col;
  int32 *order = 0;

  /*
    Order vector for col array.
  */
  order = alloc_mem( int32, nEdge );
  /*
    Initialize order. Permute columns (in old row block positions).
  */
  iord = 0;
  for (ir = 0; ir < nNod; ir++) {
    iold = perm[ir];
    for (ic = row[iold]; ic < (row[iold+1]); ic++) {
      order[ic] = iord;
      iord++;
      col[ic] = permI[col[ic]];
    }
  }
  /*
    Permute row pointers and counts.
    Use permI as tmp buffer for new row item counts.
  */
  for (ir = 0; ir < nNod; ir++) permI[ir] = (row[perm[ir]+1] - row[perm[ir]]);
  row[0] = 0;
  for (ir = 0; ir < nNod; ir++) {
    row[ir+1] = row[ir] + permI[ir];
//      printf( "%d\n", row[ir] );
  }

  if (row[nNod] != nEdge) {
    errput( "original graph was not stripped?? (%d != %d)", row[nNod], nEdge );
    ERR_CheckGo( ret );
  }
  /*
    Permute row blocks of col and val according to the order array.
  */
  for (iord = 0; iord < nEdge; iord++) {
    if (order[iord] == iord) continue;

    tmpcol = col[iord];
    neword = order[iord];

    /*
      Permutes one permutation cycle.
    */
    while (neword != iord) {
      auxint = col[neword];
      col[neword] = tmpcol;
      tmpcol = auxint;

      /*
	This gets next neword and marks old neword position in order as done.
      */
      auxint = order[neword];
      order[neword] = neword;
      neword = auxint;
    }
    /*
      Finish the cycle.
    */
    col[iord] = tmpcol;
    order[iord] = neword;
  }
  
 end_label:
  free_mem( order );

  if (ret != RET_OK) {
    errput( "graph permutation not done!" );
  }

  return( ret );
}

/*

void main(void)
{
}
*/
