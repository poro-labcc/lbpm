/**
 * 
 * Library used to calculate EDT
 * 
**/


/** operators : Basic arithmetic operation using INFTY numbers
 * 
 * David Coeurjolly (david.coeurjolly@liris.cnrs.fr) - Sept. 2004
 *
**/

#ifndef __OPERATORS_H
#define __OPERATORS_H


#define INFTY 100000001

// The sum of a and b handling INFTY
long sum(long a, long b);

//The product of a and b handling INFTY
long prod(long a, long b);

//The opposite of a  handling INFTY
long opp (long a);

// The division (integer) of divid out of divis handling INFTY
long intdivint (long divid, long divis);





////////// Functions F and Sep for the SDT labelling
/** 
 **************************************************
 * @b F
 * @param x 
 * @param i 
 * @param gi2 
 * @return Definition of a parabola
 **************************************************/
long F(int x, int i, long gi2)
{
  return sum((x-i)*(x-i), gi2);
}

/** 
 **************************************************
 * @b Sep
 * @param i 
 * @param u 
 * @param gi2 
 * @param gu2 
 * @return The absciss of the intersection point between two parabolas
 **************************************************/
long Sep(int i, int u, long gi2, long gu2) {
  return intdivint(sum( sum((long) (u*u - i*i),gu2), opp(gi2) ), 2*(u-i));
}
//////////












/////////Basic functions to handle operations with INFTY

/** 
 **************************************************
 * @b sum
 * @param a Long number with INFTY
 * @param b Long number with INFTY
 * @return The sum of a and b handling INFTY
 **************************************************/
long sum(long a, long b) 
{
  if ((a==INFTY) || (b==INFTY))     
    return INFTY;    
  else 
    return a+b;
}

/** 
 **************************************************
 * @b prod
 * @param a Long number with INFTY
 * @param b Long number with INFTY
 * @return The product of a and b handling INFTY
 **************************************************/
long prod(long a, long b) 
{
  if ((a==INFTY) || (b==INFTY)) 
    return INFTY;  
  else 
    return a*b;
}
/** 
 **************************************************
 * @b opp
 * @param a Long number with INFTY
 * @return The opposite of a  handling INFTY
 **************************************************/
long opp (long a) {
  if (a == INFTY) {
    return INFTY;
  }
  else {
    return -a;
  }
}

/** 
 **************************************************
 * @b intdivint
 * @param divid Long number with INFTY
 * @param divis Long number with INFTY
 * @return The division (integer) of divid out of divis handling INFTY
 **************************************************/
long intdivint (long divid, long divis) {
  if (divis == 0) 
    return  INFTY;
  if (divid == INFTY) 
    return  INFTY;
  else 
    return  divid / divis;
}


#endif
