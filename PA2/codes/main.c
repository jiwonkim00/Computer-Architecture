#include <stdio.h>

/*
 * CS:APP Data Lab
 *
 * <2019-11563 김지원>
 *
 * bits.c - Source file with your solutions to the Lab.
 *          This is the file you will hand in to your instructor.
 *
 * WARNING: Do not include the <stdio.h> header; it confuses the dlc
 * compiler. You can still use printf for debugging without including
 * <stdio.h>, although you might get a compiler warning. In general,
 * it's not good practice to ignore compiler warnings, but in this
 * case it's OK.
 */
#if 0
/*
 * Instructions to Students:
 *
 * STEP 1: Read the following instructions carefully.
 */

You will provide your solution to the Data Lab by
editing the collection of functions in this source file.

INTEGER CODING RULES:

  Replace the "return" statement in each function with one
  or more lines of C code that implements the function. Your code
  must conform to the following style:

  int Funct(arg1, arg2, ...) {
      /* brief description of how your implementation works */
      int var1 = Expr1;
      ...
      int varM = ExprM;

      varJ = ExprJ;
      ...
      varN = ExprN;
      return ExprR;
  }

  Each "Expr" is an expression using ONLY the following:
  1. Integer constants 0 through 255 (0xFF), inclusive. You are
      not allowed to use big constants such as 0xffffffff.
  2. Function arguments and local variables (no global variables).
  3. Unary integer operations ! ~
  4. Binary integer operations & ^ | + << >>

  Some of the problems restrict the set of allowed operators even further.
  Each "Expr" may consist of multiple operators. You are not restricted to
  one operator per line.

  You are expressly forbidden to:
  1. Use any control constructs such as if, do, while, for, switch, etc.
  2. Define or use any macros.
  3. Define any additional functions in this file.
  4. Call any functions.
  5. Use any other operations, such as &&, ||, -, or ?:
  6. Use any form of casting.
  7. Use any data type other than int.  This implies that you
     cannot use arrays, structs, or unions.


  You may assume that your machine:
  1. Uses 2s complement, 32-bit representations of integers.
  2. Performs right shifts arithmetically.
  3. Has unpredictable behavior when shifting if the shift amount
     is less than 0 or greater than 31.


EXAMPLES OF ACCEPTABLE CODING STYLE:
  /*
   * pow2plus1 - returns 2^x + 1, where 0 <= x <= 31
   */
  int pow2plus1(int x) {
     /* exploit ability of shifts to compute powers of 2 */
     return (1 << x) + 1;
  }

  /*
   * pow2plus4 - returns 2^x + 4, where 0 <= x <= 31
   */
  int pow2plus4(int x) {
     /* exploit ability of shifts to compute powers of 2 */
     int result = (1 << x);
     result += 4;
     return result;
  }

FLOATING POINT CODING RULES

For the problems that require you to implement floating-point operations,
the coding rules are less strict.  You are allowed to use looping and
conditional control.  You are allowed to use both ints and unsigneds.
You can use arbitrary integer and unsigned constants. You can use any arithmetic,
logical, or comparison operations on int or unsigned data.

You are expressly forbidden to:
  1. Define or use any macros.
  2. Define any additional functions in this file.
  3. Call any functions.
  4. Use any form of casting.
  5. Use any data type other than int or unsigned.  This means that you
     cannot use arrays, structs, or unions.
  6. Use any floating point data types, operations, or constants.


NOTES:
  1. Use the dlc (data lab checker) compiler (described in the handout) to
     check the legality of your solutions.
  2. Each function has a maximum number of operations (integer, logical,
     or comparison) that you are allowed to use for your implementation
     of the function.  The max operator count is checked by dlc.
     Note that assignment ('=') is not counted; you may use as many of
     these as you want without penalty.
  3. Use the btest test harness to check your functions for correctness.
  4. Use the BDD checker to formally verify your functions
  5. The maximum number of ops for each function is given in the
     header comment for each function. If there are any inconsistencies
     between the maximum ops in the writeup and in this file, consider
     this file the authoritative source.

/*
 * STEP 2: Modify the following functions according the coding rules.
 *
 *   IMPORTANT. TO AVOID GRADING SURPRISES:
 *   1. Use the dlc compiler to check that your solutions conform
 *      to the coding rules.
 *   2. Use the BDD checker to formally verify that your solutions produce
 *      the correct answers.
 */


#endif
/* Copyright (C) 1991-2020 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <https://www.gnu.org/licenses/>.  */
/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */
/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */
/* wchar_t uses Unicode 10.0.0.  Version 10.0 of the Unicode Standard is
   synchronized with ISO/IEC 10646:2017, fifth edition, plus
   the following additions from Amendment 1 to the fifth edition:
   - 56 emoji characters
   - 285 hentaigana
   - 3 additional Zanabazar Square characters */


/* 1
 * onebitParity - returns 1 if x contains an odd number of 0's
 *   Examples: onebitParity(5) = 0, onebitParity(7) = 1
 *   Legal ops: & ^ << >>
 *   Max ops: 20
 *   Rating: 1
 */
int onebitParity(int x) { //divide-and-conquer
    x=x^(x>>16);
    x=x^(x>>8);
    x=x^(x>>4);
    x=x^(x>>2);
    x=x^(x>>1);
    return x&1;
}


/* 2
 * checkSubstraction - Determine if can compute x-y without overflow
 *   Example: checkSubstraction(0x80000000,0x80000000) = 1,
 *            checkSubstraction(0x80000000,0x70000000) = 0,
 *   Legal ops: ! ~ & ^ | + << >>
 *   Max ops: 20
 *   Rating: 1
 */
int checkSubstraction(int x, int y) { //two cases of overflow : (+)-(-)=(-), (-)-(+)=(+)

    int x_sign = (x >> 31) & 0x1;
    int y_sign = (y >> 31) & 0x1;
    int result_sign = (x + ~y + 1) >> 31 & 0x1;

    return !((x_sign ^ y_sign) & (x_sign ^ result_sign));
    //if x, y sign is different && x, z sign is different - overflow

}

/* 3
 * twoscom2SignedVal - Convert from two's complement to signed-magnitude
 *   where the MSB is the sign bit
 *   You can assume that x > TMin
 *   Example: twoscom2SignedVal(-5) = 0x80000005.
 *   Legal ops: ! ~ & ^ | + << >>
 *   Max ops: 15
 *   Rating: 1
 */
int twoscom2SignedVal(int x) { //positive: do nothing, negative: invert all bits and add 1, after set sign bit to 1.

    int sign_bit = (x >> 31) & 0x1; // Extract the sign bit (0 for positive, 1 for negative)
    // If the sign bit is 1 (negative), negate the number and add 1
    int magnitude = (x ^ (x >> 31)) + sign_bit;
    int result = (sign_bit << 31) | (magnitude); //add the sign bit
    return result;
}


/* 4
 * nibbleReverse - Reverse nibbles(4bits) in a 32-bit word
 *   Examples: nibbleReverse(0x80000002) = 0x20000008
 *             nibbleReverse(0x89ABCDEF) = 0xFEDCBA98
 *   Legal ops: ! ~ & ^ | + << >>
 *   Max ops: 25
 *   Rating: 1
 */
int nibbleReverse(int x) {
    // Swap the nibbles using bitwise operations
    //This line extracts the least significant nibble (the rightmost 4 bits) of x and shifts it to the leftmost position (bit positions 28-31), effectively reversing its position.

    unsigned int mask1 = (0xFF | (0xFF << 8)) | ((0xFF | (0xFF << 8)) << 16);
    mask1 = mask1 & (0x0F | (0x0F << 8) | (0x0F << 16) | (0x0F << 24)); // mask1 = 0x0F0F0F0F

    unsigned int mask2 = ~mask1; // mask2 = 0xF0F0F0F0

    // Separating the nibbles
    unsigned int nibbles1 = x & mask1;
    unsigned int nibbles2 = x & mask2;

    // Reversing the nibbles and merging
    unsigned int reversed = (nibbles1 << 4) | (nibbles2 >> 4);

    return reversed;
}


/* 5
 * bitFilter - Generate a mask consisting of all 1's and filter input with it.
 * Examples: bitFilter(0xFF00, 11, 4) = 0x0F00,
 * bitFilter(0x2A00, 13, 9) = 0x2A00,
 * bitFilter(0x1300, 4, 2) = 0
 * Assume 0 <= lowbit <= 31, and 0 <= highbit <= 31
 * If lowbit > highbit, then mask should be all 0's
 *   Legal ops: & | << >>
 *   Max ops: 20
 *   Rating: 1
 */
int bitFilter(int input, int highbit, int lowbit) {
    int highbit_mask;
    int temp;
    input = (input >> lowbit) << lowbit;    //remove lowbit amount. to 0

    highbit_mask = (1 << highbit) ;
    highbit_mask |= (highbit_mask >> 1);
    highbit_mask |= (highbit_mask >> 2);
    highbit_mask |= (highbit_mask >> 4);
    highbit_mask |= (highbit_mask >> 8);
    highbit_mask |= (highbit_mask >> 16);   //make highbit mask that consists of highbit+1 number of 1's

    return (input & highbit_mask);

}
/*
 * first, move x right and left shift by lowbit - makes the low bit part all 0.
 * create a mask for highbit and do and with the x.
 * Creating a mask for high bit : do 1 << highbit, do >>1, OR these two, Do this for 2, 4, 8 until you create a mask
 */


/* 6
 * addAndDivideBy4 - adds two numbers and divide by 4 (round toward 0). But when overflow occurs
 *          while adding two numbers, returns the first operand.
 *   Examples: addAndDivideBy4(1073741824,1073741824) = 1073741824
 *             addAndDivideBy4(-2147483648,-1) = -2147483648
 *             addAndDivideBy4(32,9) = 10
 *             addAndDivideBy4(-22,9) = -3
 *   Legal ops: ! ~ & ^ | + << >>
 *   Max ops: 20
 *   Rating: 1
 */
int addAndDivideBy4(int x, int y) {
    int sum = x + y;
    int overflow = ((sum ^ x) & (sum ^ y)) >> 31;

    int sum_div_4 = sum >> 2;

    // Handle rounding towards zero for negative numbers
    int mask = (sum & (1 << 31)) >> 30;
    sum_div_4 += mask & ((sum & 3) >> 1);

    return (overflow & x) | (~overflow & sum_div_4);
    //    int sum = x + y;
//
//    // Perform the division (right shift by 2)
//    int result = sum >> 2;
//
//    // If the sum is negative, increment the result by 1
//    result = result + ((sum >> 31) & 1);
//
//    // Calculate overflow detection
//    int overflow = ((x ^ sum) & (y ^ sum)) >> 31;
//
//    // Use bitwise operations to handle overflow and return the appropriate value
//    return (overflow & x) | (~overflow & result);

    //    // Perform the addition
//    int sum = x + y;
//    // Calculate the result of dividing by 4 with rounding towards zero
//
//    int sign = (sum >> 31) & 1; // store sign bit
//    int offset = sign & 3; // If negative, add 3 to offset rounding
//    int result = (sum + offset) >> 2;
//
//    // Calculate overflow detection
//    int overflow = ((x ^ sum) & (y ^ sum)) >> 31;
//    // Use bitwise operations to handle overflow
//    return (overflow & x) | (~overflow & result);

}


/* 7
 * numZerosFirst - returns count of number of continuous 0's from first bits
 *   Example: numZerosFirst(0) = 32
 *   Example: numZerosFirst(0x80000000) = 0
 *   Example: numZerosFirst(0x40000000) = 1
 *   Example: numZerosFirst(0x00008000) = 16
 *   Legal ops: ! ~ & ^ | + << >>
 *   Max ops: 50
 *   Rating: 1
 */
int numZerosFirst(int x) {
    int result;
    int mask;
    int isZero = !x;

    /* Check the top 16 bits */
    mask = !(x >> 16) << 4; /* If top 16 bits are 0, mask = 16, else 0 */
    x = x << mask; /* Shift left by mask bits if there were zeros */
    result = mask; /* t */

    /* Check the top 8 bits of the current value */
    mask = !(x >> 24) << 3;
    x = x << mask;
    result += mask;

    /* Check the top 4 bits */
    mask = !(x >> 28) << 2;
    x = x << mask;
    result += mask;

    /* Check the top 2 bits */
    mask = !(x >> 30) << 1;
    x = x << mask;
    result += mask;

    /* Check the top bit */
    result += !(x >> 31);
    result = result + (isZero << 1) + (~isZero + 1);
    return result;

}


/* 8
 * absFloat - Return bit-level equivalent of absolute value of f for
 *   floating point argument f.
 *   Both the argument and result are passed as unsigned int's, but
 *   they are to be interpreted as the bit-level representations of
 *   single-precision floating point values.
 *   When argument is NaN, return argument..
 *   Legal ops: Any integer/unsigned operations incl. ||, &&. also if, while
 *   Max ops: 10
 *   Rating: 1
 */
unsigned absFloat(unsigned uf) {
    int zf = uf & 0x7fffffff;       // Clear the sign bit to get the absolute value
    if(zf > 0x7f800000) {   // check if NaN
        return uf;
    } else {
        return zf;
    }}


/* 9
 * castFloat2Int - Return bit-level equivalent of expression (int) f
 *   for floating point argument f.
 *   Argument is passed as unsigned int, but
 *   it is to be interpreted as the bit-level representation of a
 *   single-precision floating point value.
 *   Anything out of range (including NaN and infinity) should return
 *   0x80000000u.
 *   Legal ops: Any integer/unsigned operations incl. ||, &&. also if, while
 *   Max ops: 30
 *   Rating: 1
 */
int castFloat2Int(unsigned uf) {
    int sign = uf >> 31 & 0x1;  //extract sign bit (bit 31)
    int exp = uf >> 23 & 0xff;  //extract exponent (bits 30 to 23)
    int frac = uf & 0x7fffff;   //extracts the fraction (bits 22 to 0)

    if (exp == 0 && frac == 0) {    //checks if a zero value, return 0
        return 0;
    }
    else if (exp < 126) {   //checks underflow, returns 0 in this case
        return 0;
    }
    else if (exp > 157) {   //checks overflow, returns the minimum possible signed integer
        return 0x80000000;
    }
    else {  //uf is a normalized floating-point number, so proceed with the conversion
        exp -= 127; //calculate effective exponent by subtracting 127 from the extracted exponent.
        frac |= 0x800000; // Add implied 1 bit

        if (exp > 23) {
            exp -= 23;
            frac <<= exp;
        } else {
            exp = 23 - exp;
            frac >>= exp;
        }   //simulates reversing the normalization shift by shifting the fraction left if exp > 23 or right if exp < 23

        if (sign) {
            return ~frac + 1;
        } else {
            return frac;
        }   //If the sign bit is 1 (indicating a negative number), it takes the two's complement of the fraction to obtain the negative integer. If the sign bit is 0, it returns the fraction as an unsigned integer.
    }
}


/* 10
 * compareFloat - Compute f < g for floating point arguments f and g.
 *   Both the arguments are passed as unsigned int's, but
 *   they are to be interpreted as the bit-level representations of
 *   single-precision floating point values.
 *   If either argument is NaN, return 0.
 *   +0 and -0 are considered equal.
 *   Legal ops: Any integer/unsigned operations incl. ||, &&. also if, while
 *   Max ops: 30
 *   Rating: 1
 */
int compareFloat(unsigned uf, unsigned ug) {
    unsigned int sign_f = uf >> 31;
    unsigned int exp_f = (uf >> 23) & 0xFF;
    unsigned int frac_f = uf & 0x7FFFFF;

    unsigned int sign_g = ug >> 31;
    unsigned int exp_g = (ug >> 23) & 0xFF;
    unsigned int frac_g = ug & 0x7FFFFF;

    /* NaN check for f or g */
    if ((exp_f == 0xFF && frac_f != 0) || (exp_g == 0xFF && frac_g != 0)) {
        return 0;
    }

    /* Considering +0 and -0 as equal */
    if ((exp_f == 0 && frac_f == 0) && (exp_g == 0 && frac_g == 0)) {
        return 0;
    }

    /* Adjusting the sign bit to make comparison easier */
    if (sign_f) {
        uf = ~uf + 1;
    }

    if (sign_g) {
        ug = ~ug + 1;
    }

    return uf < ug;
}


int main() {
    int x = 0x80000000;
    int y = 0x7ffffffe;
//    int mask0 = (0x0F << 4*0) | (0x0F << 4*2) | (0x0F << 4*4) | (0x0F << 4*6);
//    int mask1 = mask0 << 4;
//    printf("%d  %d ", mask0, mask1 );


//    int highbit = 11;
//    int highbit_mask;
//    highbit_mask = (1 << highbit) ;
//    highbit_mask |= (highbit_mask >> 1);
//    highbit_mask |= (highbit_mask >> 2);
//    highbit_mask |= (highbit_mask >> 4);
//    highbit_mask |= (highbit_mask >> 8);
//    highbit_mask |= (highbit_mask >> 16);
//    printf("%d ", highbit_mask);
    int sum = x+y;
    int sign = ((x+y)>>31 & 0x1);
    int offset = ((sign<<1) | sign) & 3;
    int before = (sum + offset);
    int result = before >> 2;
    printf("sum: %d, sign: %d, offset: %d, before shift: %d, after shift: %d", sum, sign, offset, before, result);
    printf("%d", 1 & 3);
    //printf("Answer : %d", addAndDivideBy4(0x80000000, 0));
    return 0;

    /*
     * Example: numZerosFirst(0) = 32
 *   Example: numZerosFirst(0x80000000) = 0
 *   Example: numZerosFirst(0x40000000) = 1
 *   Example: numZerosFirst(0x00008000) = 16
     */
}
