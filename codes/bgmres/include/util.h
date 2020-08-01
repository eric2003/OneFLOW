#ifndef UTILROUTINE
#define UTILROUTINE


/** *********************************************************************************
 * @file util.h
 * @class ArrayUtils
 * @author Kelly Black <kjblack@gmail.com>
 * @version 0.1
 * @copyright BSD 2-Clause License
 *
 * @section LICENSE
 *
 * Copyright (c) 2014, Kelly Black
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 *
 * Class to provide a set of basic utilities that are used by numerous
 * other classes.
 *
 * This is the definition (header) file for the ArrayUtils class. It
 * includes the definitions for the methods that are used to construct
 * and delete arrays used in a variety of other classes.
 *
 *
 * @brief Header file for the basic utilities associated with managing
 * arrays.
 *
 * ********************************************************************************* */


template <class number>
class ArrayUtils
{

public:
	
/** ************************************************************************
 * Base constructor  for the ArrayUtils class.
 *
 * There is not anything to do so this is an empty method.
 *
 *  ************************************************************************ */
	ArrayUtils(){};  

	// Define the methods used to allocate the memory for and define
	// multi-dimensional arrays.
	static number *****fivetensor(int n1,int n2,int n3,int n4,int n5); 
	static number ****fourtensor(int n1,int n2,int n3,int n4);
	static number ***threetensor(int n1,int n2,int n3);
	static number **twotensor(int n1,int n2);
	static number *onetensor(int n1);

	// Define the methods used to delete the memory that was allocated
	// for multi-dimensional arrays.
	static void delfivetensor(number *****u);
	static void delfourtensor(number ****u);
	static void delthreetensor(number ***u);
	static void deltwotensor(number **u);
	static void delonetensor(number *u);

};


#include "util.cpp"


#endif

