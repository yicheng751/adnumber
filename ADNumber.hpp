#ifndef ADDNUMBER_HPP
#define	ADDNUMBER_HPP

/*!
 *  Software to compute derivatives. Support for forward, reverse, partial, 
 *  nth, and nth partial derivatives.
 */

/*!
 *   This library is dual-licensed: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version <N> as 
 *   published by the Free Software Foundation. For the terms of this 
 *   license, see licenses/gpl_v<N>.txt or <http://www.gnu.org/licenses/>.
 *
 *   You are free to use this library under the terms of the GNU General
 *   Public License, but WITHOUT ANY WARRANTY; without even the implied 
 *   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *   See the GNU General Public License for more details.
 *
 *  Alternatively, you can license this library under a commercial
 *  license.
 *
 *               ADNumber Commercial License(ACL)
 *               ================================
 * ---------------------------------------------------------------------
 * 
 *  The license is non-exclusively granted to a single person or company,
 *  per payment of the license fee, for the lifetime of that person or
 *  company. The license is non-transferrable.
 * 
 *  The ACL grants the licensee the right to use ADNumber in commercial
 *  closed-source projects. Modifications may be made to ADNumber with no
 *  obligation to release the modified source code. ADNumber (or pieces
 *  thereof) may be included in any number of projects authored (in whole
 *  or in part) by the licensee.
 * 
 *  The licensee may use any version of ADNumber, past, present or future,
 *  as is most convenient. This license does not entitle the licensee to
 *  receive any technical support, updates or bug fixes, except as such are
 *  made publicly available to all ADNumber users.
 * 
 *  The licensee may not sub-license ADNumber itself, meaning that any
 *  commercially released product containing all or parts of ADNumber must
 *  have added functionality beyond what is available in ADNumber;
 *  ADNumber itself may not be re-licensed by the licensee.
 * 
 *  To obtain a commercial license agreement, contact:
 * 
 *  Matthew Supernaw
 *  msupernaw@gmail.com
 * 
 */

/*!
 *
 * File:   ADNumber.hpp
 * Author: Matthew R. Supernaw
 *
 * Created on January 4, 2012, 7:08 PM
 */



/*!
 *Modification history
 * 
 * 2/28/2013 - fixed bug in operator -(T lhs, ADNumber rhs) M.S.
 * 
 *
 */

//#define USE_MEMORY_POOL need to create a thread safe memory pool like GCPool

#include <complex>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <limits>
#include <stdint.h>
#include <assert.h>


#define USE_CLFMALLOC//USE_MEMORY_POOL //USE_THREAD_ALLOCATOR

#ifdef USE_THREAD_ALLOCATOR

#include "threadalloc/threadalloc.h"
//#define REDEFINE_DEFAULT_NEW_OPERATOR

#endif
//#define USE_CLFMALLOC
//USE_TS_MEMORY_POOL
//USE_TCMALLOC 
/// USE_CLFMALLOC //USE_MEMORY_POOL 
//USE_CLFMALLOC

#ifdef USE_CLFMALLOC
#include "clfmalloc.h"
#endif

#ifdef USE_TS_MEMORY_POOL
#include "ts_memory_pool/static_mem_pool.h"
#endif


#ifdef USE_MEMORY_POOL
#include "memory_pool/memory_pool.hpp"
#endif


#ifdef AD_DEBUG
#define AD_TRACE(x) std::cout<<x;
#else
#define AD_TRACE(x) //std::cout<<x;
#endif




namespace ad {

    /*!
     * Creates a unique identifier.
     * @return 
     */
    class IDGenerator {
    public:
        static IDGenerator * instance();

        const uint32_t next() {
            return _id++;
        }
    private:

        IDGenerator() : _id(1) {
        }

        uint32_t _id;
    };


    static IDGenerator * only_copy;

    inline IDGenerator *
    IDGenerator::instance() {
        if (!only_copy) {
            only_copy = new IDGenerator();
        }
        return only_copy;
    }

    //!Operations used by ExpressionTree

    enum Operation {
        MINUS = 0,
        PLUS,
        MULTIPLY,
        DIVIDE,
        SIN,
        COS,
        TAN,
        ASIN,
        ACOS,
        ATAN,
        ATAN2, //atan(adnumber,adnumber)
        ATAN3, //atan(T,adnumber)
        ATAN4, //atan(adnumber,T)
        SQRT,
        POW, //pow(adnumber,adnumber)
        POW1, //pow(T,adnumber)
        POW2, //pow(adnumber,T)
        LOG,
        LOG10,
        EXP,
        SINH,
        COSH,
        TANH,
        ABS,
        FABS,
        FLOOR,
        CONSTANT,
        VARIABLE,
        NONE
    };

    /*!
     *
     * Template class Expression Tree. Used to record evaluated expression.
     * Ultimately used for arbitrary order reverse mode differentiation,
     * uncertainty, and calculation of propagated error.
     * 
     */
    template<class T>
    class Expression {
    public:

        enum MODE {
            FORWARD = 0,
            REVERSE
        };

        /*!
         * Constructor.
         * 
         * @param id
         * @param op
         * @param value
         * @param left
         * @param right
         */
        Expression(const unsigned int &id, const Operation &op,
                const T &value,
                Expression<T>* left,
                Expression<T>* right) :
        id_m(id),
        op_m(op),
        value_m(value),
        left_m(left),
        right_m(right),
        epsilon_m(std::numeric_limits<T>::epsilon()) {

        }

        /*!
         * Default constructor.
         */
        Expression() :
        left_m(NULL),
        right_m(NULL),
        value_m(T(1.0)),
        id_m(0),
        op_m(VARIABLE),
        epsilon_m(std::numeric_limits<T>::epsilon()) {

        }

        /*!
         *Copy Constructor.
         */
        Expression(const Expression &orig) :
        id_m(orig.id_m),
        op_m(orig.op_m),
        value_m(orig.value_m),
        left_m(NULL),
        right_m(NULL),
        epsilon_m(std::numeric_limits<T>::epsilon()) {

            if (orig.left_m != NULL) {
                this->left_m = orig.left_m->Clone();
            }
            if (orig.right_m != NULL) {
                this->right_m = orig.right_m->Clone();
            }

            // this = orig.Clone();
        }

        /*!
         * Destructor.
         */
        virtual ~Expression() {

            if (this->left_m) {
                delete this->left_m;
            }

            if (this->right_m) {
                delete this->right_m;
            }
        }


#ifdef USE_CLFMALLOC

        void* operator new (size_t size) throw (std::bad_alloc) {
            assert(size == sizeof (Expression));
            void* ptr = malloc(size);
            return ptr;
        }

        void operator delete (void* ptr)throw () {
            free(ptr);
        }
#endif

        size_t Size() {
            return this->Size(this);
        }

        /*!
         * Create a clone of this expression. The same as using
         * copy constructor.
         * @return 
         */
        Expression<T>* Clone() {

            return new Expression(*this);
        }

        /*!
         * Evaluate this expression. 
         * @return 
         */
        const T Evaluate(MODE mode = REVERSE) const {

            T l = T(0);
            T r = T(0);

            switch (mode) {
                case FORWARD:

                    if (this->left_m != NULL) {
                        l = this->left_m->Evaluate(mode);
                    }


                    if (this->right_m != NULL) {
                        r = this->right_m->Evaluate(mode);
                    }


                    break;
                case REVERSE:

                    if (this->right_m != NULL) {
                        r = this->right_m->Evaluate(mode);
                    }

                    if (this->left_m != NULL) {
                        l = this->left_m->Evaluate(mode);
                    }

                    break;

                default:
                    if (this->right_m != NULL) {
                        r = this->right_m->Evaluate(mode);
                    }

                    if (this->left_m != NULL) {
                        l = this->left_m->Evaluate(mode);
                    }

                    break;

            }


            switch (op_m) {
                case CONSTANT:
                    AD_TRACE("CONST[" << this->value_m << "]")
                    return this->value_m;
                case VARIABLE:
                    AD_TRACE("VAR[" << this->value_m << "]")
                    return this->value_m;
                case MINUS:
                    AD_TRACE(" - ")
                    return (l - r);
                case PLUS:
                    AD_TRACE(" + ")
                    return (l + r);
                case DIVIDE:
                    AD_TRACE(" / ")
                    return (l / r);
                case MULTIPLY:
                    AD_TRACE(" * ")
                    return (l * r);
                case SIN:
                    AD_TRACE(" sin ")
                    return sin(l);
                case COS:
                    AD_TRACE(" cos ")
                    return cos(l);
                case TAN:
                    AD_TRACE(" tan ")
                    return tan(l);
                case ASIN:
                    AD_TRACE(" asin ")
                    return asin(l);
                case ACOS:
                    AD_TRACE(" acos ")
                    return acos(l);
                case ATAN:
                    AD_TRACE(" atan ")
                    return atan(l);
                case ATAN2:
                    AD_TRACE(" atan2 ")
                    return atan2(l, r);
                    //                case ATAN3:
                    //                    break;
                    //                case ATAN4:
                    //                    break;
                case SQRT:
                    AD_TRACE(" sqrt ")
                    return sqrt(l);
                case POW:
                    AD_TRACE(" pow ")
                    return pow(l, r);
                    //                case POW1:
                    //                    break;
                    //                case POW2:
                    //                    break;
                case LOG:
                    AD_TRACE(" log ")
                    return log(l);
                case LOG10:
                    AD_TRACE(" log10 ")
                    return log10(l);
                case EXP:
                    AD_TRACE(" exp ")
                    return exp(l);
                case SINH:
                    AD_TRACE(" sinh ")
                    return sinh(l);
                case COSH:
                    AD_TRACE(" cosh ")
                    return cosh(l);
                case TANH:
                    AD_TRACE(" tanh ")
                    return tanh(l);
                case FABS:
                    AD_TRACE(" fabs ")
                    return fabs(l);
                case ABS:
                    AD_TRACE(" abs ")
                    return abs(l);
                case FLOOR:
                    AD_TRACE(" floor ")
                    return floor(l);
                case NONE:
                    AD_TRACE(" none ")
                    return this->value_m;
                default:
                    return T(0);
            }
            return T(0);
        }

        /*!
         * Evaluate the proagated error for this expression. 
         * @return 
         */
        const T PropagatedError() const {

            T a; // = this->left_->Evaluate();
            T err_a; // = this->left_->Error();
            T b; // = this->right_->Evaluate();
            T err_b; // = this->right_->Error();
            T temp_;
            T error;

            switch (op_m) {


                case CONSTANT:
                    return this->epsilon_m;
                case VARIABLE:
                    return this->epsilon_m;
                case MINUS:
                    return (this->left_m->PropagatedError() * this->left_m->PropagatedError() + this->right_m->PropagatedError() * this->right_m->PropagatedError());
                case PLUS:
                    return (this->left_m->PropagatedError() * this->left_m->PropagatedError() + this->right_m->PropagatedError() * this->right_m->PropagatedError());
                case DIVIDE:

                    if (this->left_m) {
                        a = this->left_m->Evaluate();
                    } else {
                        a = T(0);
                    }
                    err_a = this->left_m->PropagatedError();
                    if (this->right_m) {
                        b = this->right_m->Evaluate();
                    } else {
                        b = T(0);
                    }
                    err_b = this->right_m->PropagatedError();
                    return std::sqrt((err_a * err_a) / (b * b) + (a * a) *(err_b * err_b)*(b * b * b * b));

                case MULTIPLY:

                    if (this->left_m) {
                        a = this->left_m->Evaluate();
                    } else {
                        a = T(0);
                    }
                    err_a = this->left_m->PropagatedError();
                    if (this->right_m) {
                        b = this->right_m->Evaluate();
                    } else {
                        b = T(0);
                    }
                    err_b = this->right_m->PropagatedError();

                    return std::sqrt((b * b)* (err_a * err_a)+(a * a)*(err_b * err_b));

                case SIN:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();

                    return std::fabs(std::cos(a) * err_a);
                case COS:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();

                    return std::fabs(std::sin(a) * err_a);
                case TAN:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();
                    return std::fabs(err_a / std::pow(std::cos(a), T(2.0)));
                case ASIN:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();
                    return (err_a / sqrt(T(1.0) - pow(a, T(2.0))));
                case ACOS:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();
                    return (err_a / sqrt(T(1.0) + pow(a, 2.0)));
                case ATAN:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();
                    return (err_a / sqrt(T(1.0) + pow(a, 2.0)));
                case ATAN2:

                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();
                    b = this->right_m->Evaluate();
                    err_b = this->right_m->PropagatedError();

                    temp_ = fabs(a / b) *
                            sqrt(pow(err_a / a, T(2.0)) +
                            pow(err_b / b, T(2.0)));
                    error = fabs(temp_ / (T(1.0) + pow(a / b, T(2.0))));

                    return error;
                case ATAN3:
                    break;
                case ATAN4:
                case SQRT:

                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();

                    return std::fabs(err_a / T(2) * sqrt(a));
                case POW:

                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();
                    b = this->right_m->Evaluate();
                    err_b = this->right_m->PropagatedError();

                    return sqrt((b * b) * std::pow(a, (b - T(1)))*(err_a * err_a) + std::log(a) * std::log(a) * std::pow(a, b)*(err_b * err_b));
                case POW1:
                    break;
                case POW2:
                    break;
                case LOG:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();
                    return std::fabs(err_a / a);
                case LOG10:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();
                    return std::fabs(err_a / a * std::log10(T(10.0)));
                case EXP:
                    if (this->left_m) {
                        a = this->left_m->Evaluate();
                        b = this->left_m->Evaluate();
                    } else {
                        a = std::numeric_limits<T>::epsilon();
                        b = std::numeric_limits<T>::epsilon();
                    }

                    return std::fabs(std::exp(a) * b);
                case SINH:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();

                    return std::fabs(std::cosh(a) * err_a);
                case COSH:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();

                    return std::fabs(std::sinh(a) * err_a);
                case TANH:
                    a = this->left_m->Evaluate();
                    err_a = this->left_m->PropagatedError();
                    return std::fabs(err_a / std::pow(std::cosh(a), T(2.0)));
                case FABS:
                    return this->left_m->PropagatedError();
                case ABS:
                    return this->left_m->PropagatedError();
                case FLOOR:
                    return this->left_m->PropagatedError();
                case NONE:
                    return this->left_m->PropagatedError();
                default:
                    return this->epsilon_m;
            }
            return this->epsilon_m;
        }

        /*!
         * Builds a expression tree representing the derivative of this
         * tree.(reverse mode) 
         * 
         * @return Expression<T>
         */
        Expression<T>* Differentiate() {


            Expression<T>* ret = new Expression<T > ();


            switch (op_m) {

                case CONSTANT:
                    //f(x) = C
                    //f'(x) = 0

                    ret->op_m = CONSTANT;
                    ret->value_m = T(0.0);

                    //                 
                    return ret;

                case VARIABLE:
                    //f(x) = x
                    //f'(x) = 1

                    ret->op_m = CONSTANT;
                    ret->value_m = T(1.0);


                    return ret;
                case MINUS:
                    //f(x) = g(x) - h(x)
                    //f'(x) = g'(x) - h'(x)

                    ret->op_m = MINUS;
                    if (this->left_m != NULL) {
                        ret->left_m = this->left_m->Differentiate();

                    }

                    if (this->right_m != NULL) {
                        ret->right_m = this->right_m->Differentiate();
                    }

                    return ret;
                case PLUS:
                    //f(x) = g(x) + h(x)
                    //f'(x) = g'(x) + h'(x)

                    ret->op_m = PLUS;
                    if (this->left_m != NULL) {
                        ret->left_m = this->left_m->Differentiate();
                    }

                    if (this->right_m != NULL) {
                        ret->right_m = this->right_m->Differentiate();
                    }


                    return ret;
                case DIVIDE:
                    //f(x) = g(x)/h(x);
                    //f'(x) = (g'(x)h(x) - g(x)h'(x))/h(x)^2

                    ret->op_m = DIVIDE;

                    ret->left_m = new Expression<T > (); //g'(x)h(x) - g(x)h'(x)
                    ret->left_m->op_m = MINUS;


                    ret->left_m->left_m = new Expression<T > (); //g'(x)h(x)
                    ret->left_m->left_m->op_m = MULTIPLY;
                    if (this->left_m != NULL) {
                        ret->left_m->left_m->left_m = this->left_m->Differentiate();
                    }
                    ret->left_m->left_m->right_m = this->right_m->Clone();

                    ret->left_m->right_m = new Expression<T > (); //g(x)h'(x)
                    ret->left_m->right_m->op_m = MULTIPLY;
                    ret->left_m->right_m->left_m = this->left_m->Clone();
                    if (this->right_m != NULL) {
                        ret->left_m->right_m->right_m = this->right_m->Differentiate();
                    }


                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = MULTIPLY;
                    ret->right_m->left_m = this->right_m->Clone();
                    ret->right_m->right_m = this->right_m->Clone();


                    return ret;

                case MULTIPLY:
                    //f(x) = g(x)h(x);
                    //f'(x) = g'(x)h(x) + g(x)h'(x)

                    if (this->left_m->op_m == CONSTANT
                            && this->right_m->op_m != CONSTANT) {
                        ret->op_m = MULTIPLY;

                        ret->left_m = this->left_m->Clone();
                        ret->right_m = this->right_m->Differentiate();


                    } else if (this->right_m->op_m == CONSTANT
                            && this->left_m->op_m != CONSTANT) {
                        ret->op_m = MULTIPLY;

                        ret->left_m = this->left_m->Differentiate();
                        ret->right_m = this->right_m->Clone();

                    } else {



                        ret->op_m = PLUS;

                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;

                        ret->left_m->right_m = this->right_m->Clone();

                        if (this->right_m != NULL) {
                            ret->left_m->left_m = this->left_m->Differentiate();
                        }
                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = MULTIPLY;

                        ret->right_m->left_m = this->left_m->Clone();
                        if (this->left_m != NULL) {
                            ret->right_m->right_m = this->right_m->Differentiate();
                        }



                    }

                    return ret;

                case SIN:
                    //f'(x) = cos(x)

                    ret->op_m = COS;
                    ret->left_m = this->left_m->Clone();


                    return ret;

                case COS:
                    //f'(x) = -sin(x)

                    ret->op_m = MULTIPLY;


                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = CONSTANT;
                    ret->left_m->value_m = T(-1.0);

                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = SIN;
                    ret->right_m->left_m = this->left_m->Clone();


                    return ret;
                case TAN:
                    //f(x) = tan(x)
                    //f'(x) = (1/cos(x))(1/cos(x))

                    ret->op_m = MULTIPLY;

                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = DIVIDE;

                    ret->left_m->left_m = new Expression<T > ();
                    ret->left_m->left_m->op_m = CONSTANT;
                    ret->left_m->left_m->value_m = T(1.0);

                    ret->left_m->right_m = new Expression<T > ();
                    ret->left_m->right_m->op_m = COS;
                    ret->left_m->right_m->left_m = this->left_m->Clone();

                    ret->right_m = ret->left_m->Clone();


                    return ret;
                case ASIN:
                    //f(x) = asin(x)
                    //f'(x) = 1/(2 sqrt(1-x^2)= 1/(pow((1-pow(x,2)),0.5)

                    ret->op_m = DIVIDE;

                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = CONSTANT;
                    ret->left_m ->value_m = T(1.0);

                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = POW;

                    ret->right_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->op_m = MINUS;

                    ret->right_m->left_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->left_m->op_m = CONSTANT;
                    ret->right_m->left_m->left_m->value_m = T(1.0);

                    ret->right_m->left_m->right_m = new Expression<T > ();
                    ret->right_m->left_m->right_m->op_m = POW;
                    ret->right_m->left_m->right_m->left_m = this->left_m->Clone();

                    ret->right_m->left_m->right_m->right_m = new Expression<T > ();
                    ret->right_m->left_m->right_m->right_m->op_m = CONSTANT;
                    ret->right_m->left_m->right_m->right_m->value_m = T(2.0);

                    ret->right_m->right_m = new Expression<T > ();
                    ret->right_m->right_m->op_m = CONSTANT;
                    ret->right_m->right_m->value_m = T(0.5);

                    return ret;
                case ACOS:
                    //f(x) = acos(x)
                    //f'(x) = -1/(sqrt(1-x^2) = -1/(pow((1-pow(x,2)),0.5)
                    //-1/sqrt(1-x^2)

                    ret->op_m = DIVIDE;

                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = CONSTANT;
                    ret->left_m ->value_m = T(-1.0);

                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = POW;

                    ret->right_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->op_m = MINUS;

                    ret->right_m->left_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->left_m->op_m = CONSTANT;
                    ret->right_m->left_m->left_m->value_m = T(1.0);

                    ret->right_m->left_m->right_m = new Expression<T > ();
                    ret->right_m->left_m->right_m->op_m = POW;
                    ret->right_m->left_m->right_m->left_m = this->left_m->Clone();

                    ret->right_m->left_m->right_m->right_m = new Expression<T > ();
                    ret->right_m->left_m->right_m->right_m->op_m = CONSTANT;
                    ret->right_m->left_m->right_m->right_m->value_m = T(2.0);

                    ret->right_m->right_m = new Expression<T > ();
                    ret->right_m->right_m->op_m = CONSTANT;
                    ret->right_m->right_m->value_m = T(0.5);

                    return ret;
                case ATAN:
                    //f(x) = atan(x)
                    //f'(x) 1/(x^2+1)

                    ret->op_m = DIVIDE;
                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = CONSTANT;
                    ret->left_m ->value_m = T(1.0);


                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = PLUS;

                    ret->right_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->op_m = MULTIPLY;
                    ret->right_m->left_m->left_m = this->left_m->Clone();
                    ret->right_m->left_m->right_m = this->left_m->Clone();


                    ret->right_m->right_m = new Expression<T > ();
                    ret->right_m->right_m->op_m = CONSTANT;
                    ret->right_m->right_m->value_m = T(1.0);


                    return ret;
                case ATAN2:
                    //f(x) = atan2(x,y)
                    //f'(x) y/(x^2+y^2)

                    ret->op_m = DIVIDE;
                    ret->left_m = this->right_m->Clone(); //y


                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = PLUS;

                    ret->right_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->op_m = MULTIPLY;
                    ret->right_m->left_m->left_m = this->left_m->Clone();
                    ret->right_m->left_m->right_m = this->left_m->Clone();


                    ret->right_m->right_m = new Expression<T > ();
                    ret->right_m->right_m->op_m = MULTIPLY;
                    ret->right_m->right_m->left_m = this->right_m->Clone();
                    ret->right_m->right_m->right_m = this->right_m->Clone();


                    return ret;

                    //  case ATAN4:
                case SQRT:
                    //f(x) = sqrt(x)
                    //f'(x) = .5/sqrt(x)

                    ret->op_m = DIVIDE;
                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = CONSTANT;
                    ret->left_m->value_m = T(0.5);

                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = SQRT;
                    ret->right_m->left_m = this->left_m->Clone();



                    return ret;
                case POW:
                    //f(x) =  x^y
                    //f'(x) = yx^y-1

                    ret->op_m = MULTIPLY;

                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = MULTIPLY;
                    ret->left_m->left_m = this->left_m->Differentiate();
                    ret->left_m->right_m = this->right_m->Clone();



                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = POW;


                    ret->right_m->left_m = this->left_m->Clone();


                    ret->right_m->right_m = new Expression<T > ();
                    ret->right_m->right_m->op_m = MINUS;
                    ret->right_m->right_m->left_m = this->right_m->Clone();

                    ret->right_m->right_m->right_m = new Expression<T > ();
                    ret->right_m->right_m->right_m->op_m = CONSTANT;
                    ret->right_m->right_m->right_m->value_m = T(1.0);



                    return ret;
                    //                case POW1:
                    //
                    //                    break;
                    //                    //                return pow(this->left_value_, this->right_value_ - T(1.0));
                    //                case POW2:
                    //                    break;
                    //                    //                return pow(this->left_value_, this->right_value_ - T(1.0));
                case LOG:
                    //f(x) = log(x)
                    //f'(x) = 1/x

                    ret->op_m = DIVIDE;
                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = CONSTANT;
                    ret->left_m->value_m = T(1.0);

                    ret->right_m = this->left_m->Clone();



                    return ret;
                case LOG10:
                    //f(x) = log10(x)
                    //f'(x) = 1/(xlog(10))

                    ret->op_m = DIVIDE;
                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = CONSTANT;
                    ret->left_m->value_m = T(1.0);

                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = MULTIPLY;

                    ret->right_m->left_m = this->left_m->Clone();

                    ret->right_m->right_m = new Expression<T > ();
                    ret->right_m->right_m->op_m = CONSTANT;
                    ret->right_m->right_m->value_m = log(T(10.0));
                    /*
                    ret->right_->right_->left_ = new Expression<T > ();
                    ret->right_->right_->left_->op_ = CONSTANT;
                    ret->right_->right_->left_->value_ = T(10.0);
                     */

                    return ret;
                case EXP:
                    //f(x) = e^x
                    //f'(x) =e^x

                    ret->op_m = EXP;
                    ret->left_m = this->left_m->Clone();


                    return ret;
                case SINH:
                    //f(x) = sinh(x)
                    //f'(x) = cosh(x)

                    ret->op_m = COSH;
                    ret->left_m = this->left_m->Clone();


                    return ret;
                case COSH:
                    ret->op_m = SINH;
                    ret->left_m = this->left_m->Clone();


                    return ret;
                case TANH:
                    //f(x) = tanh(x)
                    //f'(x) sech^2


                    ret->op_m = MULTIPLY;

                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = DIVIDE;
                    ret->left_m->left_m = new Expression<T > ();
                    ret->left_m->left_m->op_m = CONSTANT;
                    ret->left_m->left_m->value_m = T(1.0);

                    ret->left_m->right_m = new Expression<T > ();
                    ret->left_m->right_m->op_m = COSH;
                    ret->left_m->right_m->left_m = this->left_m->Clone();

                    ret->right_m = ret->left_m->Clone();


                    return ret;

                case FABS:

                    ret->op_m = DIVIDE;
                    ret->left_m = this->left_m->Clone();

                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = FABS;
                    ret->right_m->left_m = this->left_m->Clone();



                    return ret;
                case ABS:

                    ret->op_m = DIVIDE;
                    ret->left_m = this->left_m->Clone();

                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = ABS;
                    ret->right_m->left_m = this->left_m->Clone();


                    return ret;

                case FLOOR:

                    ret->op_m = FLOOR;
                    ret->left_m = this->left_m->Clone();


                    return ret;
                case NONE://shouldn't happen.
                    return this->Clone();

                default:
                    return NULL;
            }
            return NULL;
        }

        bool HasID(const uint32_t &id) {
            //     std::cout << this->id_ << " ?= " << id << "\n";
            if (this->id_m == id) {
                return true;
            }
            if (this->left_m) {
                if (this->left_m->HasID(id)) {
                    return true;
                }
            }

            if (this->right_m) {
                if (this->right_m->HasID(id)) {
                    return true;
                }
            }

            return false;
        }

        /*!
         * Builds a expression tree representing the derivative with respect to 
         * some ADNumber via its id.(reverse mode) 
         * 
         * @return Expression<T>
         */
        Expression<T>* Differentiate(const uint32_t &id) {
            //#warning need to check partial derivatives....

            Expression<T>* ret = new Expression<T > ();


            switch (op_m) {

                case CONSTANT:
                    //f(x) = C
                    //f'(x) = 0

                    ret->op_m = CONSTANT;
                    ret->value_m = T(0); //this->value_;


                    return ret;

                case VARIABLE:
                    if (this->id_m == id) {
                        //f(x) = x
                        //f'(x) = 1

                        ret->op_m = CONSTANT;
                        ret->value_m = T(1.0);


                        return ret;
                    } else {//constant
                        //f(x) = C
                        //f'(x) = 0
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);
                        return ret;
                    }
                case MINUS:

                    //f(x) = g(x) - h(x)
                    //f'(x) = g'(x) - h'(x)

                    ret->op_m = MINUS;
                    if (this->left_m) {
                        ret->left_m = this->left_m->Differentiate(id);

                    }

                    if (this->right_m) {
                        ret->right_m = this->right_m->Differentiate(id);
                    }

                    return ret;

                case PLUS:

                    //f(x) = g(x) + h(x)
                    //f'(x) = g'(x) + h'(x)

                    ret->op_m = PLUS;
                    if (this->left_m) {
                        ret->left_m = this->left_m->Differentiate(id);
                    }

                    if (this->right_m) {
                        ret->right_m = this->right_m->Differentiate(id);
                    }


                    return ret;

                case DIVIDE:

                    //f(x) = g(x)/h(x);
                    //f'(x) = (g'(x)h(x) - g(x)h'(x))/h(x)^2


                    ret->op_m = DIVIDE;

                    ret->left_m = new Expression<T > (); //g'(x)h(x) - g(x)h'(x)
                    ret->left_m->op_m = MINUS;


                    ret->left_m->left_m = new Expression<T > (); //g'(x)h(x)
                    ret->left_m->left_m->op_m = MULTIPLY;
                    if (this->left_m) {
                        ret->left_m->left_m->left_m = this->left_m->Differentiate(id);
                    }
                    ret->left_m->left_m->right_m = this->right_m->Clone();

                    ret->left_m->right_m = new Expression<T > (); //g(x)h'(x)
                    ret->left_m->right_m->op_m = MULTIPLY;
                    ret->left_m->right_m->left_m = this->left_m->Clone();
                    if (this->right_m) {
                        ret->left_m->right_m->right_m = this->right_m->Differentiate(id);
                    }


                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = MULTIPLY;
                    ret->right_m->left_m = this->right_m->Clone();
                    ret->right_m->right_m = this->right_m->Clone();


                    return ret;

                case MULTIPLY:
                    //f(x) = g(x)h(x);
                    //f'(x) = g'(x)h(x) + g(x)h'(x)

                    if (this->left_m->op_m == CONSTANT
                            && this->right_m->op_m != CONSTANT) {
                        ret->op_m = MULTIPLY;
                        if (this->left_m) {
                            ret->left_m = this->left_m->Clone();
                        }
                        if (this->right_m) {
                            ret->right_m = this->right_m->Differentiate(id);
                        }


                    } else if (this->right_m->op_m == CONSTANT
                            && this->left_m->op_m != CONSTANT) {
                        ret->op_m = MULTIPLY;
                        if (this->left_m) {
                            ret->left_m = this->left_m->Differentiate(id);
                        }
                        if (this->right_m) {
                            ret->right_m = this->right_m->Clone();
                        }
                    } else {



                        ret->op_m = PLUS;

                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;

                        ret->left_m->right_m = this->right_m->Clone();

                        if (this->right_m != NULL) {
                            ret->left_m->left_m = this->left_m->Differentiate(id);
                        }
                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = MULTIPLY;

                        ret->right_m->left_m = this->left_m->Clone();
                        if (this->left_m != NULL) {
                            ret->right_m->right_m = this->right_m->Differentiate(id);
                        }



                    }
                    return ret;

                case SIN:

                    if (this->left_m->HasID(id)) {
                        //f'(x) = cos(x)

                        ret->op_m = MULTIPLY;
                        ret->left_m = this->left_m->Differentiate(id);
                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = COS;
                        ret->right_m->left_m = this->left_m->Clone();

                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }

                case COS:
                    if (this->left_m->HasID(id)) {
                        //f'(x) = -sin(x)

                        ret->op_m = MULTIPLY;


                        ret->left_m = this->left_m->Differentiate(id);
                        ret->right_m = new Expression<T > ();

                        ret->right_m->op_m = MULTIPLY;
                        ret->right_m->left_m = new Expression<T > ();
                        ret->right_m->left_m->op_m = CONSTANT;
                        ret->right_m->left_m->value_m = T(-1.0);

                        ret->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->op_m = SIN;
                        ret->right_m->right_m->left_m = this->left_m->Clone();


                        return ret;

                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case TAN:
                    if (this->left_m->HasID(id)) {
                        //f'(x) = 1/cos(x)

                        ret->op_m = MULTIPLY;
                        ret->left_m = this->left_m->Differentiate(id);


                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = MULTIPLY;

                        ret->right_m->left_m = new Expression<T > ();
                        ret->right_m->left_m->op_m = DIVIDE;


                        ret->right_m->left_m->left_m = new Expression<T > ();
                        ret->right_m->left_m->left_m->op_m = CONSTANT;
                        ret->right_m->left_m->left_m->value_m = T(1.0);


                        ret->right_m->left_m->right_m = new Expression<T > ();
                        ret->right_m->left_m->right_m->op_m = COS;
                        ret->right_m->left_m->right_m->left_m = this->left_m->Clone();


                        ret->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->op_m = DIVIDE;


                        ret->right_m->right_m->left_m = new Expression<T > ();
                        ret->right_m->right_m->left_m->op_m = CONSTANT;
                        ret->right_m->right_m->left_m->value_m = T(1.0);


                        ret->right_m->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->right_m->op_m = COS;
                        ret->right_m->right_m->right_m->left_m = this->left_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case ASIN:

                    if (this->left_m->HasID(id)) {
                        //f(x) = asin(x)
                        //f'(x) = 1/(2 sqrt(1-x^2)= 1/(pow((1-pow(x,2)),0.5)

                        ret->op_m = MULTIPLY;
                        ret->left_m = this->left_m->Differentiate(id);


                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = DIVIDE;

                        ret->right_m->left_m = new Expression<T > ();
                        ret->right_m->left_m->op_m = CONSTANT;
                        ret->right_m->left_m ->value_m = T(1.0);

                        ret->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->op_m = POW;

                        ret->right_m->right_m->left_m = new Expression<T > ();
                        ret->right_m->right_m->left_m->op_m = MINUS;

                        ret->right_m->right_m->left_m->left_m = new Expression<T > ();
                        ret->right_m->right_m->left_m->left_m->op_m = CONSTANT;
                        ret->right_m->right_m->left_m->left_m->value_m = T(1.0);

                        ret->right_m->right_m->left_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->left_m->right_m->op_m = POW;
                        ret->right_m->right_m->left_m->right_m->left_m = this->left_m->Clone();

                        ret->right_m->right_m->left_m->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->left_m->right_m->right_m->op_m = CONSTANT;
                        ret->right_m->right_m->left_m->right_m->right_m->value_m = T(2.0);

                        ret->right_m->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->right_m->op_m = CONSTANT;
                        ret->right_m->right_m->right_m->value_m = T(0.5);

                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case ACOS:

                    if (this->left_m->HasID(id)) {
                        //f(x) = acos(x)
                        //f'(x) = -1/(sqrt(1-x^2) = -1/(pow((1-pow(x,2)),0.5)
                        //-1/sqrt(1-x^2)
                        ret->op_m = MULTIPLY;
                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;
                        ret->left_m->left_m = new Expression<T > ();

                        ret->left_m->left_m->op_m = CONSTANT;
                        ret->left_m->left_m->value_m = T(-1.0);


                        ret->left_m->right_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = DIVIDE;

                        ret->right_m->left_m = new Expression<T > ();
                        ret->right_m->left_m->op_m = CONSTANT;
                        ret->right_m->left_m ->value_m = T(1.0);

                        ret->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->op_m = POW;

                        ret->right_m->right_m->left_m = new Expression<T > ();
                        ret->right_m->right_m->left_m->op_m = MINUS;

                        ret->right_m->right_m->left_m->left_m = new Expression<T > ();
                        ret->right_m->right_m->left_m->left_m->op_m = CONSTANT;
                        ret->right_m->right_m->left_m->left_m->value_m = T(1.0);

                        ret->right_m->right_m->left_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->left_m->right_m->op_m = POW;
                        ret->right_m->right_m->left_m->right_m->left_m = this->left_m->Clone();

                        ret->right_m->right_m->left_m->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->left_m->right_m->right_m->op_m = CONSTANT;
                        ret->right_m->right_m->left_m->right_m->right_m->value_m = T(2.0);

                        ret->right_m->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->right_m->op_m = CONSTANT;
                        ret->right_m->right_m->right_m->value_m = T(0.5);

                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case ATAN:
                    if (this->left_m->HasID(id)) {
                        //f(x) = atan(x)
                        //f'(x) 1/(x^2+1)

                        ret->op_m = DIVIDE;
                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;
                        ret->left_m->right_m = new Expression<T > ();

                        ret->left_m->right_m->op_m = CONSTANT;
                        ret->left_m->right_m->value_m = T(1.0);


                        ret->left_m->left_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = PLUS;

                        ret->right_m->left_m = new Expression<T > ();
                        ret->right_m->left_m->op_m = MULTIPLY;
                        ret->right_m->left_m->left_m = this->left_m->Clone();
                        ret->right_m->left_m->right_m = this->left_m->Clone();


                        ret->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->op_m = CONSTANT;
                        ret->right_m->right_m->value_m = T(1.0);


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;

                    }
                case ATAN2:
                    //if w.r.t. check both expressions for id
                    if (this->left_m->HasID(id)) {
                        //f(x) = atan2(x,y)
                        //f'(x) y/(x^2+y^2)

                        ret->op_m = DIVIDE;
                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;
                        ret->left_m->left_m = this->right_m->Clone(); //y
                        ret->left_m->right_m = left_m->Differentiate(id);


                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = PLUS;

                        ret->right_m->left_m = new Expression<T > ();
                        ret->right_m->left_m->op_m = MULTIPLY;
                        ret->right_m->left_m->left_m = this->left_m->Clone();
                        ret->right_m->left_m->right_m = this->left_m->Clone();


                        ret->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->op_m = MULTIPLY;
                        ret->right_m->right_m->left_m = this->right_m->Clone();
                        ret->right_m->right_m->right_m = this->right_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case ATAN3:

                    //can be removed.
                    break;

                case ATAN4:
                    break;
                case SQRT:
                    if (this->left_m->HasID(id)) {
                        //f(x) = sqrt(x)
                        //f'(x) = .5/sqrt(x)

                        ret->op_m = DIVIDE;
                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;

                        ret->left_m->right_m = new Expression<T > ();
                        ret->left_m->right_m->value_m = T(0.5);

                        ret->left_m->left_m = this->left_m->Differentiate(id);


                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = SQRT;
                        ret->right_m->left_m = this->left_m->Clone();

                        //std::cout<<ret->ToString();

                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case POW:

                    if (this->left_m->HasID(id)) {
                        //f(x) =  x^y
                        //f'(x) = yx^y-1

                        ret->op_m = MULTIPLY;

                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;
                        ret->left_m->left_m = this->left_m->Differentiate(id);
                        ret->left_m->right_m = this->right_m->Clone();



                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = POW;


                        ret->right_m->left_m = this->left_m->Clone();


                        ret->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->op_m = MINUS;
                        ret->right_m->right_m->left_m = this->right_m->Clone();

                        ret->right_m->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->right_m->op_m = CONSTANT;
                        ret->right_m->right_m->right_m->value_m = T(1.0);



                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                    //                case POW1:
                    //
                    //                    break;
                    //                    //                return pow(this->left_value_, this->right_value_ - T(1.0));
                    //                case POW2:
                    //                    break;
                    //                    //                return pow(this->left_value_, this->right_value_ - T(1.0));
                case LOG:
                    if (this->left_m->HasID(id)) {
                        //f(x) = log(x)
                        //f'(x) = 1/x

                        ret->op_m = DIVIDE;
                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;
                        ret->left_m->left_m = new Expression<T > ();
                        ret->left_m->left_m->op_m = CONSTANT;
                        ret->left_m->left_m->value_m = T(1.0);
                        ret->left_m->right_m = this->left_m->Differentiate(id);

                        ret->right_m = this->left_m->Clone();



                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case LOG10:
                    //f(x) = log10(x)
                    //f'(x) = 1/(xlog(10))

                    if (this->left_m->HasID(id)) {



                        ret->op_m = DIVIDE;

                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;

                        ret->left_m->left_m = new Expression<T > ();
                        ret->left_m->left_m->op_m = CONSTANT;
                        ret->left_m->left_m->value_m = T(1.0);

                        ret->left_m->right_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = MULTIPLY;

                        ret->right_m->left_m = this->left_m->Clone();

                        ret->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->op_m = CONSTANT;
                        ret->right_m->right_m->value_m = log(T(10.0));


                        return ret;
                    } else {
                        ret->op_m = LOG;
                        ret->left_m = this->Clone();


                        return ret;
                    }
                case EXP:
                    //f(x) = e^x
                    //f'(x) =e^x

                    if (this->left_m->HasID(id)) {

                        ret->op_m = MULTIPLY;
                        ret->left_m = this->left_m->Differentiate(id);


                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = EXP;
                        ret->right_m->left_m = this->left_m->Clone();



                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case SINH:
                    if (this->left_m->HasID(id)) {
                        //f(x) = sinh(x)
                        //f'(x) = cosh(x)

                        ret->op_m = MULTIPLY;
                        ret->left_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = COSH;
                        ret->right_m->left_m = this->left_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case COSH:
                    if (this->left_m->HasID(id)) {

                        ret->op_m = MULTIPLY;
                        ret->left_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = SINH;
                        ret->right_m->left_m = this->left_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case TANH:
                    //f(x) = tanh(x)
                    //f'(x) =1- tanh(x)*tanh(x)


                    if (this->left_m->HasID(id)) {

                        ret->op_m = MULTIPLY;

                        ret->left_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = MULTIPLY;
                        ret->right_m->left_m = new Expression<T > ();


                        ret->right_m->left_m->op_m = DIVIDE;
                        ret->right_m->left_m->left_m = new Expression<T > ();
                        ret->right_m->left_m->left_m->op_m = CONSTANT;
                        ret->right_m->left_m->left_m->value_m = T(1.0);


                        ret->right_m->left_m->right_m = new Expression<T > ();
                        ret->right_m->left_m->right_m->op_m = COSH;
                        ret->right_m->left_m->right_m->left_m = this->left_m->Clone();


                        ret->right_m->right_m = ret->right_m->left_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }

                case FABS:

                    if (this->left_m->HasID(id)) {

                        ret->op_m = DIVIDE;
                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;

                        ret->left_m->left_m = this->left_m->Differentiate(id);
                        ret->left_m->right_m = this->left_m->Clone();


                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = FABS;
                        ret->right_m->left_m = this->left_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case FLOOR:
                    if (this->left_m->id_m == id) {



                        ret->op_m = MULTIPLY;

                        ret->left_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = FLOOR;
                        ret->right_m->left_m = this->left_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case NONE://shouldn't happen.
                    return this->Clone();

                default:
                    return NULL;
            }
            return NULL;
        }


#define SPEED_UP

#ifdef SPEED_UP

        /**
         * Returns the evaluated derivative of this expression tree. While the
         * derivative is computed, no expression tree manipulations are made.
         * @param id
         * @return 
         */
        T EvaluateDerivative(const uint32_t &id) {
            //#warning need to check partial derivatives....

            T ret, g, h = T(-999.0);



            switch (op_m) {

                case CONSTANT:
                    //f(x) = C
                    //f'(x) = 0

                    return T(0);

                case VARIABLE:
                    if (this->id_m == id) {
                        //f(x) = x
                        //f'(x) = 1


                        return T(1.0);
                    } else {//constant
                        //f(x) = C
                        //f'(x) = 0

                        return T(0.0);
                    }
                case MINUS:

                    //f(x) = g(x) - h(x)
                    //f'(x) = g'(x) - h'(x)


                    return this->left_m->EvaluateDerivative(id) - this->right_m->EvaluateDerivative(id);

                case PLUS:

                    //f(x) = g(x) + h(x)
                    //f'(x) = g'(x) + h'(x)


                    return this->left_m->EvaluateDerivative(id) + this->right_m->EvaluateDerivative(id);


                case DIVIDE:

                    //f(x) = g(x)/h(x);
                    //f'(x) = (g'(x)h(x) - g(x)h'(x))/h(x)^2


                    ret = (this->left_m->EvaluateDerivative(id) * this->right_m->Evaluate() -
                            this->left_m->Evaluate() * this->right_m->EvaluateDerivative(id)) /
                            (this->right_m->Evaluate() * this->right_m->Evaluate());


                    return ret;

                case MULTIPLY:
                    //f(x) = g(x)h(x);
                    //f'(x) = g'(x)h(x) + g(x)h'(x)

                    if (this->left_m->op_m == CONSTANT
                            && this->right_m->op_m != CONSTANT) {

                        ret = this->left_m->Evaluate() * this->right_m->EvaluateDerivative(id);


                    } else if (this->right_m->op_m == CONSTANT
                            && this->left_m->op_m != CONSTANT) {

                        ret = this->left_m->EvaluateDerivative(id) * this->right_m->Evaluate();
                    } else {

                        //g'(x)h(x) + g(x)h'(x)

                        ret = this->left_m->EvaluateDerivative(id) * this->right_m->Evaluate() +
                                this->left_m->Evaluate() * this->right_m->EvaluateDerivative(id);


                    }
                    return ret;

                case SIN:

                    if (this->left_m->HasID(id)) {
                        //f'(x) = cos(x)
                        ret = this->left_m->EvaluateDerivative(id) *
                                std::cos(this->left_m->Evaluate());

                        return ret;
                    } else {
                        return T(0.0);
                    }

                case COS:
                    if (this->left_m->HasID(id)) {
                        //f'(x) = -sin(x)


                        g = this->left_m->EvaluateDerivative(id);

                        ret = g * T(-1.0) * std::sin(this->left_m->Evaluate());

                        return ret;

                    } else {

                        return T(0.0);
                    }
                case TAN:
                    if (this->left_m->HasID(id)) {
                        //f'(x) = 1/cos(x)


                        g = this->left_m->EvaluateDerivative(id);

                        ret = g * ((T(1.0) / std::cos(this->left_m->Evaluate()))*(T(1.0) / std::cos(this->left_m->Evaluate())));


                        return ret;
                    } else {

                        return T(0.0);
                    }
                case ASIN:

                    if (this->left_m->HasID(id)) {


                        //f(x) = asin(x)
                        //f'(x) = 1/(2 sqrt(1-x^2)= 1/(pow((1-pow(x,2)),0.5)


                        g = this->left_m->EvaluateDerivative(id);

                        ret = (g * T(1.0) / std::pow((T(1.0) - std::pow(this->left_m->Evaluate(), T(2.0))), T(0.5)));

                        return ret;
                    } else {
                        return T(0.0);
                    }
                case ACOS:

                    if (this->left_m->HasID(id)) {
                        g = this->left_m->EvaluateDerivative(id);

                        ret = (g * T(-1.0) / std::pow((T(1.0) - std::pow(this->left_m->Evaluate(), T(2.0))), T(0.5)));

                        return ret;
                    } else {

                        return T(0.0);
                    }
                case ATAN:
                    if (this->left_m->HasID(id)) {
                        g = this->left_m->EvaluateDerivative(id);
                        ret = (g * T(1.0) / (this->left_m->Evaluate() * this->left_m->Evaluate() + T(1.0)));

                        return ret;
                    } else {
                        //                        ret->op_m = CONSTANT;
                        //                        ret->value_m = T(0.0);
                        return T(0.0);

                    }
                case ATAN2:
                    //if w.r.t. check both expressions for id
                    if (this->left_m->HasID(id)) {
                        //f(x) = atan2(x,y)
                        //f'(x) y/(x^2+y^2)

                        g = this->left_m->EvaluateDerivative(id);
                        ret = (this->right_m->Evaluate() * g / (this->left_m->Evaluate() * this->left_m->Evaluate()+(this->right_m->Evaluate() * this->right_m->Evaluate())));

                        return ret;
                    } else {

                        return T(0.0);
                    }
                case ATAN3:

                    //can be removed.
                    break;

                case ATAN4:
                    break;
                case SQRT:
                    if (this->left_m->HasID(id)) {
                        //f(x) = sqrt(x)
                        //f'(x) = .5/sqrt(x)
                        g = this->left_m->EvaluateDerivative(id);
                        ret = g * T(.5) / std::sqrt(this->left_m->Evaluate());


                        return ret;
                    } else {
                        return T(0.0);
                    }
                case POW:

                    if (this->left_m->HasID(id)) {
                        //f(x) =  x^y
                        //f'(x) = yx^y-1
                        ret = (this->left_m->EvaluateDerivative(id) * this->right_m->Evaluate()) *
                                std::pow(this->left_m->Evaluate(), (this->right_m->Evaluate() - T(1.0)));

                        return ret;
                    } else {

                        return T(0.0);
                    }

                case LOG:
                    if (this->left_m->HasID(id)) {
                        //f(x) = log(x)
                        //f'(x) = 1/x
                        ret = (this->left_m->EvaluateDerivative(id) * T(1.0)) / this->left_m->Evaluate();
                   
                        return ret;
                    } else {
                                              
                        return T(0.0);
                    }
                case LOG10:
                    //f(x) = log10(x)
                    //f'(x) = 1/(xlog(10))

                    if (this->left_m->HasID(id)) {

                        ret = (this->left_m->EvaluateDerivative(id) * T(1.0)) / (this->left_m->Evaluate() * std::log(T(10.0)));

                        return ret;
                    } else {
                         return T(0.0);
                    }
                case EXP:
                    //f(x) = e^x
                    //f'(x) =e^x

                    if (this->left_m->HasID(id)) {
                        ret = this->left_m->EvaluateDerivative(id) * std::exp(this->left_m->Evaluate());
                     
                        return ret;
                    } else {
                        
                        return T(0.0);
                    }
                case SINH:
                    if (this->left_m->HasID(id)) {
                        //f(x) = sinh(x)
                        //f'(x) = cosh(x)
                        return this->left_m->EvaluateDerivative(id) * std::cosh(this->left_m->Evaluate());
                       
                        return ret;
                    } else {
                       
                        return T(0.0);
                    }
                case COSH:
                    if (this->left_m->HasID(id)) {
                        return this->left_m->EvaluateDerivative(id) * std::sinh(this->left_m->Evaluate());
                      
                        return ret;
                    } else {
                       
                        return ret;
                    }
                case TANH:
                    //f(x) = tanh(x)
                    //f'(x) =1- tanh(x)*tanh(x)


                    if (this->left_m->HasID(id)) {

                        ret = this->left_m->EvaluateDerivative(id)*(T(1.0) / std::cosh(this->left_m->Evaluate()))*(T(1.0) / std::cosh(this->left_m->Evaluate()));


                        return ret;
                    } else {
                       
                        return T(0.0);
                    }

                case FABS:

                    if (this->left_m->HasID(id)) {

                        ret= (this->left_m->EvaluateDerivative(id) * this->left_m->Evaluate()) /
                                std::fabs(this->left_m->Evaluate());
                       
                        return ret;
                    } else {
                       
                        return T(0.0);
                    }
                case FLOOR:
                    if (this->left_m->id_m == id) {

                        ret= this->left_m->EvaluateDerivative(id) * std::floor(this->left_m->Evaluate());
                       
                        return ret;
                    } else {
                        
                        return ret;
                    }
                case NONE://shouldn't happen.
                    return ret;

                default:
                    return ret;
            }
            return NULL;
        }

        T FindValue(const uint32_t &id) {

            if (this->id_m == id) {
                return this->value_m;
                ;
            }
            if (this->left_m) {
                if (this->left_m->HasID(id)) {
                    return this->left_m->FindValue(id);
                }
            }

            if (this->right_m) {
                if (this->right_m->HasID(id)) {
                    return this->left_m->FindValue(id);
                }
            }

            return T(0);
        }
#endif

        /*!
         * Builds a expression tree representing the integral with respect to 
         * some ADNumber via its id.(reverse mode) 
         * 
         * @return Expression<T>
         */
        Expression<T>* Integral(const uint32_t &id) {
            //#warning need to check partial derivatives....

            Expression<T>* ret = new Expression<T > ();


            switch (op_m) {

                case CONSTANT:


                    ret->op_m = MULTIPLY;
                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = CONSTANT;
                    ret->left_m->value_m = this->value_m;

                    ret->right_m = new Expression<T > ();
                    ret->right_m->id_m = id;
                    ret->right_m->value_m = this->FindValue(id);


                    return ret;

                case VARIABLE:
                    if (this->id_m == id) {


                        ret->op_m = DIVIDE;
                        ret->left_m = new Expression<T > ();
                        ret->right_m = new Expression<T > ();

                        ret->left_m->op_m = MULTIPLY;
                        ret->left_m->left_m = this->Clone();
                        ret->left_m->right_m = this->Clone();

                        ret->right_m->op_m = CONSTANT;
                        ret->right_m->value_m = T(2);


                        return ret;
                    } else {//constant
                        //f(x) = C
                        //f'(x) = 0
                        ret->op_m = MULTIPLY;
                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = VARIABLE;
                        ret->left_m->id_m = id;
                        ret->left_m->value_m = this->FindValue(id);

                        ret->right_m = this->Clone();

                        return ret;
                    }
                case MINUS:

                    //f(x) = g(x) - h(x)
                    //f'(x) = g'(x) - h'(x)

                    ret->op_m = MINUS;
                    if (this->left_m) {
                        ret->left_m = this->left_m->Integral(id);

                    }

                    if (this->right_m) {
                        ret->right_m = this->right_m->Integral(id);
                    }

                    return ret;

                case PLUS:

                    //f(x) = g(x) + h(x)
                    //f'(x) = g'(x) + h'(x)

                    ret->op_m = PLUS;
                    if (this->left_m) {
                        ret->left_m = this->left_m->Integral(id);
                    }

                    if (this->right_m) {
                        ret->right_m = this->right_m->Integral(id);
                    }


                    return ret;

                case DIVIDE:

                    ret->op_m = DIVIDE;
                    if (this->left_m) {
                        ret->left_m = this->left_m->Integral(id);
                    }

                    if (this->right_m) {
                        ret->right_m = this->right_m->Integral(id);
                    }


                    return ret;

                case MULTIPLY:
                    ret->op_m = MULTIPLY;
                    if (this->left_m) {
                        ret->left_m = this->left_m->Integral(id);
                    }

                    if (this->right_m) {
                        ret->right_m = this->right_m->Integral(id);
                    }


                    return ret;

                case SIN:

                    ret->op_m = MULTIPLY;
                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = CONSTANT;
                    ret->left_m->value_m = T(-1);

                    ret->right_m = this->Clone();
                    ret->right_m->op_m = COS;

                    return ret;



                case COS:

                    ret = this->Clone();
                    ret->op_m = SIN;
                    return ret;

                case TAN:

                    ret->op_m = LOG;
                    ret->left_m = new Expression<T > ();
                    ret->left_m->left_m = new Expression<T > ();
                    ret->left_m->left_m->op_m = CONSTANT;
                    ret->left_m->left_m->value_m = T(1);

                    ret->left_m->right_m = this->Clone();
                    ret->left_m->right_m->op_m = COS;

                    return ret;

                case ASIN:

                    ret->op_m = PLUS;

                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = MULTIPLY;
                    ret->left_m->left_m = this->left_m->Clone();
                    ret->left_m->right_m = new Expression<T > ();
                    ret->left_m->right_m->op_m = ASIN;
                    ret->left_m->right_m->left_m = this->left_m->Clone();


                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = SQRT;
                    ret->right_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->op_m = PLUS;
                    ret->right_m->left_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->left_m->op_m = CONSTANT;
                    ret->right_m->left_m->left_m->value_m = T(1);


                    ret->right_m->left_m->right_m = new Expression<T > ();
                    ret->right_m->left_m->right_m->op_m = MULTIPLY;
                    ret->right_m->left_m->right_m->left_m = this->left_m->Clone();
                    ret->right_m->left_m->right_m->right_m = this->left_m->Clone();



                    return ret;

                case ACOS:

                    ret->op_m = MINUS;

                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = MULTIPLY;
                    ret->left_m->left_m = this->left_m->Clone();
                    ret->left_m->right_m = new Expression<T > ();
                    ret->left_m->right_m->op_m = ACOS;
                    ret->left_m->right_m->left_m = this->left_m->Clone();


                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = SQRT;
                    ret->right_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->op_m = MINUS;
                    ret->right_m->left_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->left_m->op_m = CONSTANT;
                    ret->right_m->left_m->left_m->value_m = T(1);


                    ret->right_m->left_m->right_m = new Expression<T > ();
                    ret->right_m->left_m->right_m->op_m = MULTIPLY;
                    ret->right_m->left_m->right_m->left_m = this->left_m->Clone();
                    ret->right_m->left_m->right_m->right_m = this->left_m->Clone();



                    return ret;
                case ATAN:
                    ret->op_m = MINUS;

                    ret->left_m = new Expression<T > ();
                    ret->left_m->op_m = MULTIPLY;
                    ret->left_m->left_m = this->left_m->Clone();
                    ret->left_m->right_m = new Expression<T > ();
                    ret->left_m->right_m->op_m = ATAN;
                    ret->left_m->right_m->left_m = this->left_m->Clone();


                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = MULTIPLY;
                    ret->right_m->left_m = new Expression<T > ();
                    ret->right_m->left_m->value_m = T(.5);
                    ret->right_m->right_m = new Expression<T > ();
                    ret->right_m->right_m->op_m = LOG;
                    ret->right_m->right_m->left_m = new Expression<T > ();
                    ret->right_m->right_m->left_m->op_m = PLUS;
                    ret->right_m->right_m->left_m->left_m = new Expression<T > ();
                    ret->right_m->right_m->left_m->left_m->op_m = CONSTANT;
                    ret->right_m->right_m->left_m->left_m->value_m = T(1);
                    ret->right_m->right_m->left_m->right_m = new Expression<T > ();
                    ret->right_m->right_m->left_m->right_m->op_m = MULTIPLY;
                    ret->right_m->right_m->left_m->right_m->left_m = this->left_m->Clone();
                    ret->right_m->right_m->left_m->right_m->right_m = this->left_m->Clone();

                    return ret;

                case ATAN2:

                    std::cout << "Error, still haven't implemented integral of atan2!\n";
                    exit(0);
                case ATAN3:

                    //can be removed.
                    break;

                case ATAN4:
                    break;
                case SQRT:

                    ret->op_m = MULTIPLY;
                    ret->left_m = new Expression<T > ();

                    ret->left_m->op_m = DIVIDE;
                    ret->left_m->left_m = new Expression<T > ();
                    ret->left_m->left_m->op_m = CONSTANT;
                    ret->left_m->left_m->value_m = T(2);
                    ret->left_m->right_m = new Expression<T > ();
                    ret->left_m->right_m->op_m = CONSTANT;
                    ret->left_m->right_m->value_m = T(3);

                    ret->right_m = new Expression<T > ();
                    ret->right_m->op_m = POW;
                    ret->right_m->left_m = this->Clone();
                    ret->right_m->right_m = ret->left_m->Clone();

                    return ret;

                case POW:


                    ret->op_m = DIVIDE;
                    ret->left_m = new Expression<T > ();

                    ret->left_m->op_m = POW;
                    ret->left_m->left_m = this->left_m->Clone();
                    ret->left_m->right_m = new Expression<T > ();

                    ret->left_m->right_m->op_m = PLUS;
                    ret->left_m->right_m->left_m = this->right_m->Clone();
                    ret->left_m->right_m->right_m = new Expression<T > ();
                    ret->left_m->right_m->right_m->op_m = CONSTANT;
                    ret->left_m->right_m->right_m->value_m = T(1);


                    ret->right_m = ret->left_m->right_m->Clone();


                    return ret;
                    //                case POW1:
                    //
                    //                    break;
                    //                    //                return pow(this->left_value_, this->right_value_ - T(1.0));
                    //                case POW2:
                    //                    break;
                    //                    //                return pow(this->left_value_, this->right_value_ - T(1.0));
                case LOG:
                    if (this->left_m->HasID(id)) {
                        //f(x) = log(x)
                        //f'(x) = 1/x

                        ret->op_m = DIVIDE;
                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;
                        ret->left_m->left_m = new Expression<T > ();
                        ret->left_m->left_m->op_m = CONSTANT;
                        ret->left_m->left_m->value_m = T(1.0);
                        ret->left_m->right_m = this->left_m->Differentiate(id);

                        ret->right_m = this->left_m->Clone();



                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case LOG10:
                    //f(x) = log10(x)
                    //f'(x) = 1/(xlog(10))

                    if (this->left_m->HasID(id)) {



                        ret->op_m = DIVIDE;

                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;

                        ret->left_m->left_m = new Expression<T > ();
                        ret->left_m->left_m->op_m = CONSTANT;
                        ret->left_m->left_m->value_m = T(1.0);

                        ret->left_m->right_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = MULTIPLY;

                        ret->right_m->left_m = this->left_m->Clone();

                        ret->right_m->right_m = new Expression<T > ();
                        ret->right_m->right_m->op_m = CONSTANT;
                        ret->right_m->right_m->value_m = log(T(10.0));
                       
                        return ret;
                    } else {
                        ret->op_m = LOG;
                        ret->left_m = this->Clone();


                        return ret;
                    }
                case EXP:
                    //f(x) = e^x
                    //f'(x) =e^x

                    if (this->left_m->HasID(id)) {

                        ret->op_m = MULTIPLY;
                        ret->left_m = this->left_m->Differentiate(id);


                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = EXP;
                        ret->right_m->left_m = this->left_m->Clone();



                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case SINH:
                    if (this->left_m->HasID(id)) {
                        //f(x) = sinh(x)
                        //f'(x) = cosh(x)

                        ret->op_m = MULTIPLY;
                        ret->left_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = COSH;
                        ret->right_m->left_m = this->left_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case COSH:
                    if (this->left_m->HasID(id)) {

                        ret->op_m = MULTIPLY;
                        ret->left_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = SINH;
                        ret->right_m->left_m = this->left_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case TANH:
                    //f(x) = tanh(x)
                    //f'(x) =1- tanh(x)*tanh(x)


                    if (this->left_m->HasID(id)) {

                        ret->op_m = MULTIPLY;

                        ret->left_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = MULTIPLY;
                        ret->right_m->left_m = new Expression<T > ();


                        ret->right_m->left_m->op_m = DIVIDE;
                        ret->right_m->left_m->left_m = new Expression<T > ();
                        ret->right_m->left_m->left_m->op_m = CONSTANT;
                        ret->right_m->left_m->left_m->value_m = T(1.0);


                        ret->right_m->left_m->right_m = new Expression<T > ();
                        ret->right_m->left_m->right_m->op_m = COSH;
                        ret->right_m->left_m->right_m->left_m = this->left_m->Clone();


                        ret->right_m->right_m = ret->right_m->left_m->Clone();
                       

                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }

                case FABS:

                    if (this->left_m->HasID(id)) {

                        ret->op_m = DIVIDE;
                        ret->left_m = new Expression<T > ();
                        ret->left_m->op_m = MULTIPLY;

                        ret->left_m->left_m = this->left_m->Differentiate(id);
                        ret->left_m->right_m = this->left_m->Clone();


                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = FABS;
                        ret->right_m->left_m = this->left_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case FLOOR:
                    if (this->left_m->id_m == id) {



                        ret->op_m = MULTIPLY;

                        ret->left_m = this->left_m->Differentiate(id);

                        ret->right_m = new Expression<T > ();
                        ret->right_m->op_m = FLOOR;
                        ret->right_m->left_m = this->left_m->Clone();


                        return ret;
                    } else {
                        ret->op_m = CONSTANT;
                        ret->value_m = T(0.0);


                        return ret;
                    }
                case NONE://shouldn't happen.
                    return this->Clone();

                default:
                    return NULL;
            }
            return NULL;
        }

        /*!
         * Solves an expression with respect to some 
         * ADNumber via its id.
         *
         */
        void Solve(const uint32_t &id) {

        }

     
        /*!
         * Return a list of differentiable ids in this
         * expression.
         * 
         * @param vars
         */
        void VariableIds(std::vector< uint32_t> &vars) {

            if (this->left_m != NULL) {
                this->left_m->VariableIds(vars);
            }


            if (this->right_m != NULL) {
                this->right_m->VariableIds(vars);
            }

            if (this->op_m == VARIABLE) {
                bool exists = false;
                for (size_t i = 0; i < vars.size(); i++) {
                    if (vars.at(i) == this->id_m) {
                        exists = true;
                    }
                }

                if (!exists) {
                    vars.push_back(this->id_m);
                }
            }

        }

        /*!
         * Represent this expression as a string. ADNumbers are represented
         * in wkt format by value and id. Constants are represented by value.
         *
         */
        std::string ToString(std::vector<std::string> &vars) {
            std::stringstream ss;
            std::stringstream temp;

            std::string l, r;
            bool exists = false;
            std::stringstream temps;

            if (this->left_m != NULL) {
                l = this->left_m->ToString(vars);
            }
            temp.str("");

            if (this->right_m != NULL) {
                r = this->right_m->ToString(vars);
            }
            ss << "(";

            switch (this->op_m) {
                case CONSTANT:
                    ss << this->value_m << "";
                    break;
                case VARIABLE:
                    ss << "x" << this->id_m;


                    temps << "x" << this->id_m << "/*" << this->value_m << "*/";

                    for (size_t i = 0; i < vars.size(); i++) {
                        if (vars.at(i) == temps.str()) {
                            exists = true;
                        }
                    }

                    if (!exists) {
                        vars.push_back(temps.str());
                    }

                    break;
                case MINUS:
                    ss << l << " - " << r;
                    break;
                case PLUS:
                    ss << l << " + " << r;
                    break;
                case DIVIDE:
                    ss << l << " / " << r;
                    break;
                case MULTIPLY:
                    ss << l << " * " << r;
                    break;
                case SIN:
                    ss << "sin(" << l << ")";
                    break;
                case COS:
                    ss << "cos(" << l << ")";
                    break;
                case TAN:
                    ss << "tan(" << l << ")";
                    break;
                case ASIN:
                    ss << "asin(" << l << ")";
                    break;
                case ACOS:
                    ss << "acos(" << l << ")";
                    break;
                case ATAN:
                    ss << "atan(" << l << ")";
                    break;
                case ATAN2:
                    ss << "atan(" << l << "," << r << ")";
                    break;
                case ATAN3:
                    ss << "atan(" << l << "," << r << ")";
                    break;
                case ATAN4:
                    ss << "atan(" << l << "," << r << ")";
                    break;
                case SQRT:
                    ss << "sqrt(" << l << ")";
                    break;
                case POW:
                    ss << "pow(" << l << "," << r << ")";
                    break;
                case POW1:
                    ss << "pow(" << l << "," << r << ")";
                    break;
                case POW2:
                    ss << "pow(" << l << "," << r << ")";
                    break;
                case LOG:
                    ss << "log(" << l << ")";
                    break;
                case LOG10:
                    ss << "log10(" << l << ")";
                    break;
                case EXP:
                    ss << "exp(" << l << ")";
                    break;
                case SINH:
                    ss << "sinh(" << l << ")";
                    break;
                case COSH:
                    ss << "cosh(" << l << ")";
                    break;
                case TANH:
                    ss << "tanh(" << l << ")";
                    break;
                case ABS:
                    ss << "abs(" << l << ")";
                case NONE:
                    break;
                default:
                    break;
            }

            ss << ")";

            return ss.str();

        }

      

        /*!
         * Represent this expression as a string. ADNumbers are represented
         * in wkt format by value and id. Constants are represented by value.
         *
         */
        std::string ToString(MODE mode = REVERSE) {
            std::stringstream ss;
            std::stringstream temp;

            std::string l, r;

            switch (mode) {
                case FORWARD:

                    if (this->left_m != NULL) {
                        l = this->left_m->ToString(mode);
                    }


                    if (this->right_m != NULL) {
                        r = this->right_m->ToString(mode);
                    }


                    break;
                case REVERSE:

                    if (this->right_m != NULL) {
                        r = this->right_m->ToString(mode);
                    }

                    if (this->left_m != NULL) {
                        l = this->left_m->ToString(mode);
                    }

                    break;

                default:
                    if (this->right_m != NULL) {
                        r = this->right_m->ToString(mode);
                    }

                    if (this->left_m != NULL) {
                        l = this->left_m->ToString(mode);
                    }

                    break;

            }
            //            if (this->left_m != NULL) {
            //                l = this->left_m->ToString();
            //            }
            temp.str("");

            //            if (this->right_m != NULL) {
            //                r = this->right_m->ToString();
            //            }
            ss << "(";

            switch (this->op_m) {
                case CONSTANT:
                    ss << "CONST[" << this->value_m << "]";
                    break;
                case VARIABLE:
                    ss << "VAR[" << this->value_m << ",ID[" << this->id_m << "]" << "]";
                    break;
                case MINUS:
                    ss << l << " - " << r;
                    break;
                case PLUS:
                    ss << l << " + " << r;
                    break;
                case DIVIDE:
                    ss << l << " / " << r;
                    break;
                case MULTIPLY:
                    ss << l << " * " << r;
                    break;
                case SIN:
                    ss << "sin(" << l << ")";
                    break;
                case COS:
                    ss << "cos(" << l << ")";
                    break;
                case TAN:
                    ss << "tan(" << l << ")";
                    break;
                case ASIN:
                    ss << "asin(" << l << ")";
                    break;
                case ACOS:
                    ss << "acos(" << l << ")";
                    break;
                case ATAN:
                    ss << "atan(" << l << ")";
                    break;
                case ATAN2:
                    ss << "atan(" << l << "," << r << ")";
                    break;
                case ATAN3:
                    ss << "atan(" << l << "," << r << ")";
                    break;
                case ATAN4:
                    ss << "atan(" << l << "," << r << ")";
                    break;
                case SQRT:
                    ss << "sqrt(" << l << ")";
                    break;
                case POW:
                    ss << "pow(" << l << "," << r << ")";
                    break;
                case POW1:
                    ss << "pow(" << l << "," << r << ")";
                    break;
                case POW2:
                    ss << "pow(" << l << "," << r << ")";
                    break;
                case LOG:
                    ss << "log(" << l << ")";
                    break;
                case LOG10:
                    ss << "log10(" << l << ")";
                    break;
                case EXP:
                    ss << "exp(" << l << ")";
                    break;
                case SINH:
                    ss << "sinh(" << l << ")";
                    break;
                case COSH:
                    ss << "cosh(" << l << ")";
                    break;
                case TANH:
                    ss << "tanh(" << l << ")";
                    break;
                case FABS:
                    ss << "fabs(" << l << ")";
                case FLOOR:
                    ss << "floor(" << l << ")";
                case NONE:
                    break;
                default:
                    break;
            }

            ss << ")";

            return ss.str();

        }

        /*!
         * Returns the unique identifier associated with this 
         * expression.
         */
        unsigned long GetId() const {

            return id_m;
        }

        /*!
         * Set the unique identifier associated with this expression.
         * @param id
         */
        void SetId(unsigned long id) {

            this->id_m = id;
        }

        /*!
         * Returns the left branch of this expression tree.
         * @return 
         */
        Expression<T>* GetLeft() const {

            return left_m;
        }

        /*!
         * Sets the left branch of this expression tree.
         *
         */
        void SetLeft(Expression<T> *left) {

            this->left_m = left;
        }

        /*!
         * Return the operator for this expression.
         * @return 
         */
        Operation GetOp() const {

            return op_m;
        }

        /*!
         * Set the operator for this expression.
         * 
         */
        void SetOp(const Operation &op) {

            this->op_m = op;
        }

        /*!
         * Returns the right branch of this expression tree.
         * @return 
         */
        Expression<T>* GetRight() const {

            return right_m;
        }

        /*!
         * Sets the right branch of this expression tree.
         * 
         */
        void SetRight(Expression<T> *right) {

            this->right_m = right;
        }

        /*!
         * Returns the value assigned to this 
         * expression.
         * @return 
         */
        T GetValue() const {

            return value_m;
        }

        /*!
         * Sets the value assigned to this 
         * expression.
         * 
         */
        void SetValue(T value) {
            this->value_m = value;
        }



    private:
        T epsilon_m;
        Expression<T>* left_m; //left branch
        Expression<T>* right_m; //right branch
        T value_m; //raw value
        unsigned int id_m; //unique id
        Operation op_m; //operation

        size_t Size(Expression* exp) {
            size_t count = 0;
            if (exp != NULL) {
                count = 1 + Size(exp->left_m) + Size(exp->right_m);
            }
            return count;
        }



    };

    /*!
     * class ADNumber. 
     * 
     * @brief
     * A template class to perform Automatic differentiation.
     * Supports forward and reverse mode traversal of the chain rule.
     * Supports higher order derivatives, as well as  partial and nth
     * partial derivatives with respect to a ADNumber. Overrides most functions
     * in cmath.h (all in the namespace std). Works by storing evaluated 
     * expressions in a expression tree. The expression tree can than be 
     * manipulated to give nth and nth partial derivatives. Template parameter
     * should be of floating point type, either native or arbitrary precision.
     * If arbitrary precision type is used, cmath.h functions must be 
     * overridden. 
     * 
     * @author Matthew Supernaw
     * 
     * @contact msupernaw@gmail.com
     * 
     * @date January 4, 2012
     * 
     * 
     */
    template<class T>
    class ADNumber {
    public:

        /*!
         * Default constructor.
         */
        ADNumber() :
        expression_m(new Expression<T>()),
        value_m(T(0.0)),
        fderivative_m(T(1)),
        variableName_(std::string("na")),
        id_m(IDGenerator::instance()->next()) {
            //  Lock l(this->mutex_m);
            this->Initialize();

        }


        //
        //    ADNumber(T value, T derivative = T(1)) : value_(value),
        //    fderivative_(derivative),
        //    variableName_(std::string("x")),
        //    id_(IDGenerator::instance()->next()),
        //    expression_(new Expression<T>()) {
        //        this->Initialize();
        //    }

        /*!
         * Constructor
         * @param value
         * @param derivative
         */
        ADNumber(const T &value, const T &derivative = T(1.0)) :
        expression_m(new Expression<T>()),
        value_m(value),
        fderivative_m(derivative),
        variableName_(std::string("x")),
        id_m(IDGenerator::instance()->next()) {
            //  Lock l(this->mutex_m);
            this->Initialize();

        }

        /*!
         * Constructor
         * 
         * @param name
         * @param value
         * @param derivative
         */
        ADNumber(const std::string &name, const T &value, const T &derivative = T(1)) :
        value_m(value),
        fderivative_m(derivative),
        variableName_(name),
        expression_m(new Expression<T>()),
        id_m(IDGenerator::instance()->next()) {
            //  Lock l(this->mutex_m);
            this->Initialize();
        }

        /*!
         * Copy Constructor.
         * 
         * @param orig
         */
        ADNumber(const ADNumber& orig) :
        value_m(T(0.0)),
        fderivative_m(T(1)),
        variableName_(std::string("na")),
        id_m(orig.id_m) {
            // Lock l(this->mutex_m);
            //  Lock ll(orig.mutex_m);

            this->variableName_ = orig.variableName_;
            this->value_m = orig.value_m;
            this->fderivative_m = orig.fderivative_m;
            this->id_m = orig.id_m;

            if (orig.expression_m != NULL) {

                this->expression_m = orig.expression_m->Clone();
            }
            // this->Initialize();
        }


        /*!
         * Destructor.
         */
        virtual ~ADNumber() {
            //  Lock l(this->mutex_m);

            if (this->expression_m != NULL) {
                delete this->expression_m;
            }
        }

        //    /*!
        //     * returns this value.
        //     */

        operator T() const {
            return this->value_m;
        }

        /*!
         * Returns this value.
         */
        operator ADNumber<T>() const {
            //  Lock l(this->mutex_m);
            return this;
        }

        /**
         * Return the size of the underlying expression tree.
         * 
         * @return expression size
         */
        size_t Size() {
            return this->expression_m->Size();
        }

        /**
         * In member assignment operator to set this 
         * equal to another ADNumber val.
         * 
         * @param val
         * @return ADNumber
         */
        ADNumber<T> & operator =(const ADNumber<T> &val) {
            //  Lock l(this->mutex_m);
            if (val.expression_m != NULL) {
                delete this->expression_m;
            }

            this->expression_m = val.expression_m->Clone();
            this->value_m = val.GetValue();
            this->fderivative_m = val.Forward();
            this->id_m = val.id_m;
            this->variableName_ = val.variableName_;

            return *this;
        }


        /**
         * In member assignment operator to set this value 
         * equal to val with derivative set to 1.
         * 
         * @param val
         * @return ADNumber
         */
        ADNumber<T> & operator =(const T & val) {
            //  Lock l(this->mutex_m);
            this->value_m = val;
            this->fderivative_m = T(1.0);
            this->id_m = uint32_t(IDGenerator::instance()->next());
            delete this->expression_m;
            this->expression_m = new Expression<T > ();


            this->Initialize();

            return *this;
        }

        /**
         * In member addition operator.
         * Returns ADNumber<T> with this value + rhs value &
         * this derivative + rhs derivitive.
         * 
         * @param rhs
         * @return ADNumber
         */
        const ADNumber<T> operator +(const ADNumber<T>& rhs) const {
            //  Lock l(this->mutex_m);
            ADNumber<T > ret(T(this->value_m + rhs.GetValue()),
                    this->fderivative_m + rhs.Forward());

            ret.expression_m->SetOp(PLUS);
            ret.expression_m->SetLeft(this->expression_m->Clone());
            ret.expression_m->SetRight(rhs.expression_m->Clone());

            return ret;
        }

        /**
         * In member addition operator.
         * Returns ADNumber<T> with this value + rhs value &
         * this derivative + 0.
         * @param rhs
         * 
         * @return ADNumber
         */
        const ADNumber<T> operator +(const T & rhs) const {
            //  Lock l(this->mutex_m);
            ADNumber<T > ret(T(this->value_m + rhs),
                    T(this->fderivative_m));

            ret.expression_m->SetOp(PLUS);
            ret.expression_m->SetLeft(this->expression_m->Clone());

            Expression<T> *temp = new Expression<T > ();
            temp->SetOp(CONSTANT);
            temp->SetValue(rhs);
            ret.expression_m->SetRight(temp);

            return ret;
        }

        /**
         * In member subtraction operator.
         * Returns ADNumber<T> with this value - rhs value &
         * this derivative - rhs derivitive.
         * 
         * @param rhs
         * @return ADNmuber
         */
        const ADNumber<T> operator -(const ADNumber<T>& rhs) const {
            //  Lock l(this->mutex_m);
            ADNumber<T > ret(T(this->value_m - rhs.GetValue()),
                    T(this->fderivative_m - rhs.Forward()));


            ret.expression_m->SetOp(MINUS);
            ret.expression_m->SetLeft(this->expression_m->Clone());
            ret.expression_m->SetRight(rhs.expression_m->Clone());

            return ret;
        }

        /**
         * In member subtraction operator.
         * Returns ADNumber<T>(this value - rhs value, this derivative - 0).
         * 
         * @param rhs
         * @return ADNumber
         */
        const ADNumber<T> operator -(const T & rhs) const {
            //  Lock l(this->mutex_m);
            ADNumber<T > ret(T(this->value_m - rhs),
                    T(this->fderivative_m));

            ret.expression_m->SetOp(MINUS);
            ret.expression_m->SetLeft(this->expression_m->Clone());

            Expression<T> *temp = new Expression<T > ();
            temp->SetOp(CONSTANT);
            temp->SetValue(rhs);
            ret.expression_m->SetRight(temp);

            return ret;
        }

        /*!
         * In member multiplication operator.
         * Returns ADNumber<T>(this value * rhs value,
         * this value_ * rhs derivative  + rhs value * this derivative).
         */
        const ADNumber<T> operator *(const ADNumber<T>& rhs) const {
            //  Lock l(this->mutex_m);
            ADNumber<T > ret(T(this->value_m * rhs.GetValue()),
                    T(this->value_m * rhs.Forward() +
                    rhs.GetValue() * this->fderivative_m));

            ret.expression_m->SetOp(MULTIPLY);
            ret.expression_m->SetLeft(this->expression_m->Clone());
            ret.expression_m->SetRight(rhs.expression_m->Clone());

            return ret;
        }

        /*!
         * In member multiplication operator.
         * Returns ADNumber<T>(this value * rhs value,
         * this value_ * 0  + rhs value * this derivative).
         */
        const ADNumber<T> operator *(const T & rhs) const {
            //  Lock l(this->mutex_m);
            ADNumber<T > ret(T(this->value_m * rhs),
                    T(this->value_m * T(0) + rhs * this->fderivative_m));

            ret.expression_m->SetOp(MULTIPLY);
            ret.expression_m->SetLeft(this->expression_m->Clone());

            Expression<T> *temp = new Expression<T > ();
            temp->SetOp(CONSTANT);
            temp->SetValue(rhs);
            ret.expression_m->SetRight(temp);

            return ret;
        }

        /**
         * 
         * @param rhs
         * @return 
         */
        const ADNumber<T> operator /(const ADNumber<T>& rhs) const {
            //  Lock l(this->mutex_m);
            ADNumber<T > ret(T(this->value_m / rhs.value_m),
                    T((rhs.GetValue() * this->fderivative_m -
                    this->value_m * rhs.Forward()) / (rhs.GetValue() * rhs.GetValue())));

            ret.expression_m->SetOp(DIVIDE);
            ret.expression_m->SetLeft(this->expression_m->Clone());
            ret.expression_m->SetRight(rhs.expression_m->Clone());

            return ret;
        }

        /*!
         * In member division operator.
         * Returns ADNumber<T>(this value * rhs value,
         * (this value_ * rhs derivative  - rhs value * 0)/(rhs value * rhs value)).
         */
        const ADNumber<T> operator /(const T & rhs) const {
            //  Lock l(this->mutex_m);
            ADNumber<T > ret(T(this->value_m / rhs),
                    T((rhs * this->value_m - this->value_m * 0) / (rhs * rhs)));

            ret.expression_m->SetOp(DIVIDE);
            ret.expression_m->SetLeft(this->expression_m->Clone());

            Expression<T> *temp = new Expression<T > ();
            temp->SetOp(CONSTANT);
            temp->SetValue(rhs);
            ret.expression_m->SetRight(temp);

            return ret;
        }

        /*!
         * In member addition assignment operator.
         * @param rhs
         * @return 
         */
        ADNumber<T>& operator +=(const ADNumber<T>& rhs) {

            Expression<T>* exp = new Expression<T>;
            exp->SetLeft(this->expression_m);
            exp->SetRight(rhs.expression_m->Clone());
            exp->SetValue(this->GetValue() + rhs.GetValue());
            exp->SetOp(PLUS);
            this->expression_m = exp;
            this->value_m += rhs.GetValue();
            this->fderivative_m = T(this->value_m * rhs.Forward() +
                    rhs.GetValue() * this->fderivative_m);

            return *this;
        }

        /*!
         * In member addition subtraction operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T>& operator -=(const ADNumber<T>& rhs) {
            Expression<T>* exp = new Expression<T>;
            exp->SetLeft(this->expression_m);
            exp->SetRight(rhs.expression_m->Clone());
            exp->SetValue(this->GetValue() + rhs.GetValue());
            exp->SetOp(MINUS);
            this->expression_m = exp;
            this->value_m += rhs.GetValue();
            this->fderivative_m = T(this->fderivative_m - rhs.Forward());

            return *this;
        }

        /*!
         * In member multiplication assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T>& operator *=(const ADNumber<T>& rhs) {
            //            //  Lock l(this->mutex_m);
            //            ADNumber<T> temp = *this;
            //            ADNumber<T> ret = (temp * rhs);
            //            if (ret.expression_m != NULL) {
            //                delete this->expression_m;
            //                this->expression_m = ret.expression_m->Clone();
            //            }
            //            this->value_m = ret.GetValue();
            //            this->fderivative_m = ret.Forward();
            //
            //            return ret;
            *this = *this*rhs;
            return *this;
        }

        /*!
         * In member division assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T>& operator /=(const ADNumber<T>&rhs) {
            //            //  Lock l(this->mutex_m);
            //            ADNumber<T> temp = *this;
            //            ADNumber<T> ret = (temp / rhs);
            //            if (ret.expression_m != NULL) {
            //                delete this->expression_m;
            //                this->expression_m = ret.expression_m->Clone();
            //            }
            //            this->value_m = ret.GetValue();
            //            this->fderivative_m = ret.Forward();
            //
            //            return ret;
            *this = *this*rhs;
            return *this;
        }

        /*!
         * In member addition assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T>& operator +=(const T & rhs) {
            //  Lock l(this->mutex_m);
            //            ADNumber<T> temp = *this;
            //            ADNumber<T> ret = (temp + rhs);
            //            if (ret.expression_m != NULL) {
            //                delete this->expression_m;
            //                this->expression_m = ret.expression_m->Clone();
            //            }
            //            this->value_m = ret.GetValue();
            //            this->fderivative_m = ret.Forward();
            //
            //            return ret;
            *this = *this+rhs;
            return *this;
        }

        /*!
         * In member subtraction assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T>& operator -=(const T & rhs) {
            //            //  Lock l(this->mutex_m);
            //            ADNumber<T> temp = *this;
            //            ADNumber<T> ret = (temp - rhs);
            //            if (ret.expression_m != NULL) {
            //                delete this->expression_m;
            //                this->expression_m = ret.expression_m->Clone();
            //            }
            //            this->value_m = ret.GetValue();
            //            this->fderivative_m = ret.Forward();
            //
            //            return ret;
            *this = *this-rhs;
            return *this;
        }

        /*!
         * In member multiplication assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T>& operator *=(const T & rhs) {
            //            //  Lock l(this->mutex_m);
            //            ADNumber<T> temp = *this;
            //            ADNumber<T> ret = (temp * rhs);
            //            if (ret.expression_m != NULL) {
            //                delete this->expression_m;
            //                this->expression_m = ret.expression_m->Clone();
            //            }
            //            this->value_m = ret.GetValue();
            //            this->fderivative_m = ret.Forward();
            //
            //            return ret;
            *this = *this*rhs;
            return *this;
        }

        /*!
         * In member division assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T>* operator /=(const T & rhs) {
            //            //  Lock l(this->mutex_m);
            //            ADNumber<T> temp = *this;
            //            ADNumber<T> ret = (temp / rhs);
            //            if (ret.expression_m != NULL) {
            //                delete this->expression_m;
            //                this->expression_m = ret.expression_m->Clone();
            //            }
            //            this->value_m = ret.GetValue();
            //            this->fderivative_m = ret.Forward();
            //
            //            return ret;
            *this = *this / rhs;
            return *this;
        }

        /*!
         * In member suffix increment operator.
         * 
         * @return 
         */
        ADNumber<T>& operator ++() {
            //            //  Lock l(this->mutex_m);
            //            ADNumber<T> temp = *this;
            //            ADNumber<T> ret(*this+T(1.0));
            //            this->value_m = ret.GetValue();
            //            this->fderivative_m = ret.Forward();
            //            delete this->expression_m;
            //            this->expression_m = ret.expression_m->Clone();
            //
            //            return *this;
            *this = *this+T(1.0);
            return *this;
        }

        /*!
         * In member suffix decrement operator.
         * 
         * @return 
         */
        ADNumber<T> operator --() {
            //  Lock l(this->mutex_m);
            ADNumber<T> temp = *this;
            ADNumber<T> ret(temp - T(1.0));
            this->value_m = ret.GetValue();
            this->fderivative_m = ret.Forward();
            delete this->expression_m;
            this->expression_m = ret.expression_m->Clone();

            return *this;
        }

        /*!
         * In member prefix increment operator.
         * 
         * @param 
         */
        ADNumber<T> operator ++(int) {
            //  Lock l(this->mutex_m);
            ADNumber<T> temp = *this;

            ADNumber<T> ret(temp + T(1.0));
            this->value_m = ret.GetValue();
            this->fderivative_m = ret.Forward();
            delete this->expression_m;
            this->expression_m = ret.expression_m->Clone();

            return temp;

        }

        /*!
         * In member prefix decrement operator.
         * 
         * @param 
         */
        ADNumber<T> operator --(int) {
            //  Lock l(this->mutex_m);
            ADNumber<T> temp = *this;

            ADNumber<T> ret(*this-T(1.0));
            this->value_m = ret.GetValue();
            this->fderivative_m = ret.Forward();
            delete this->expression_m;
            this->expression_m = ret.expression_m->Clone();

            return *temp;
        }

        void Reset(T val = T(0), T deriv = T(1)) {
            if (this->expression_m != NULL) {
                delete this->expression_m;
            }

            this->expression_m = new Expression<T > ();
            this->value_m = val;
            this->fderivative_m = deriv;
            this->Initialize();
        }

        /*!
         * Returns the computed value.
         */
        const T GetValue() const {
            //  //  Lock l(this->mutex_m);
            return this->value_m;
        }

        /*!
         * Returns the result of all derivatives in the expression chain via forward
         * accumulation. Forward mode automatic differentiation is accomplished by
         * the dual numbers method.
         * Source:http://en.wikipedia.org/wiki/Automatic_differentiation#Automatic_differentiation_using_dual_numbers
         */
        const T Forward() const {
            ////  Lock l(this->mutex_m);
            return this->fderivative_m;
        }

        /*!
         * Returns the result of all derivatives in the expression chain via reverse
         * accumulation.  Reverse accumulation traverses the chain rule from left
         * to right, or in the case of the computational graph,
         * from top to bottom.
         * Source: http://en.wikipedia.org/wiki/Automatic_differentiation#Reverse_accumulation
         */
        const T Reverse() const {
            //Lock l(this->mutex_m);
            if (this->expression_m == NULL) {
                return T(0);
            }

            Expression<T> * exp = this->expression_m->Differentiate();
            T ret = exp->Evaluate();
            delete exp;

            return ret;

        }

        /*!
         * Derivative with respect to var0....var1 in order.
         * @param var0
         * @return numerical derivative
         */
        const T WRT(const ADNumber<T> &var0) {
            //  Lock l(this->mutex_m);
            // Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            T ret = exp->Evaluate();
            delete exp;
            return ret;
        }

        /*!
         * Derivative with respect to var0....var2 in order.
         * @param var0
         * @param var1
         * @return numerical derivative
         */
        const T WRT(const ADNumber<T> &var0, const ADNumber<T> &var1) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());

            T ret = temp->Evaluate();
            delete temp;
            delete exp;


            return ret;
        }

        /*!
         * Derivative with respect to var0....var3 in order.
         * @param var0
         * @param var1
         * @param var2
         * @return numerical derivative
         */
        const T WRT(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            T ret = exp->Evaluate();
            delete exp;

            return ret;
        }

        /*!
         * Derivative with respect to var0....var4 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @return numerical derivative
         */
        const T WRT(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            T ret = temp->Evaluate();
            delete temp;

            return ret;
        }

        /*!
         * Derivative with respect to var0....var5 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @return numerical derivative
         */
        const T WRT(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            T ret = exp->Evaluate();
            delete exp;

            return ret;
        }

        /*!
         * Derivative with respect to var0....var6 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @param var5
         * @return numerical derivative
         */
        const T WRT(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4,
                const ADNumber<T> &var5) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            temp = exp->Differentiate(var5.GetID());
            delete exp;
            T ret = temp->Evaluate();
            delete temp;

            return ret;
        }

        /*!
         * Derivative with respect to var0....var7 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @param var5
         * @param var6
         * @return numerical derivative
         */
        const T WRT(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4,
                const ADNumber<T> &var5, const ADNumber<T> &var6) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            temp = exp->Differentiate(var5.GetID());
            delete exp;
            exp = temp->Differentiate(var6.GetID());
            delete temp;
            T ret = exp->Evaluate();
            delete exp;

            return ret;
        }

        /*!
         * Derivative with respect to var0....var8 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @param var5
         * @param var6
         * @param var7
         * @return numerical derivative
         */
        const T WRT(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4,
                const ADNumber<T> &var5, const ADNumber<T> &var6, const ADNumber<T> &var7) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            temp = exp->Differentiate(var5.GetID());
            delete exp;
            exp = temp->Differentiate(var6.GetID());
            delete temp;
            temp = exp->Differentiate(var7.GetID());
            delete exp;
            T ret = temp->Evaluate();
            delete temp;

            return ret;
        }

        /*!
         * Derivative with respect to var0....var9 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @param var5
         * @param var6
         * @param var7
         * @param var8
         * @return numerical derivative
         */
        const T WRT(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4,
                const ADNumber<T> &var5, const ADNumber<T> &var6, const ADNumber<T> &var7,
                const ADNumber<T> &var8) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            temp = exp->Differentiate(var5.GetID());
            delete exp;
            exp = temp->Differentiate(var6.GetID());
            delete temp;
            temp = exp->Differentiate(var7.GetID());
            delete exp;
            exp = temp->Differentiate(var8.GetID());
            delete temp;
            T ret = exp->Evaluate();
            delete exp;

            return ret;
        }

        /*!
         * Derivative with respect to var0....var10 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @param var5
         * @param var6
         * @param var7
         * @param var8
         * @param var9
         * @return numerical derivative
         */
        const T WRT(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4,
                const ADNumber<T> &var5, const ADNumber<T> &var6, const ADNumber<T> &var7,
                const ADNumber<T> &var8, const ADNumber<T> &var9) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            temp = exp->Differentiate(var5.GetID());
            delete exp;
            exp = temp->Differentiate(var6.GetID());
            delete temp;
            temp = exp->Differentiate(var7.GetID());
            delete exp;
            exp = temp->Differentiate(var8.GetID());
            delete temp;
            temp = exp->Differentiate(var9.GetID());
            delete exp;
            T ret = temp->Evaluate();
            delete temp;

            return ret;
        }

        /*!
         * Derivative with respect to a vector of ADNumbers 
         * {var0,var1...varn} in order.
         * @param vars
         * @return numerical derivative
         */
        const T WRT(const std::vector<ADNumber<T>*> &vars) {
            //  Lock l(this->mutex_m);
            if (this->expression_m == NULL) {
                return T(0);
            } else {
                if (vars.size() == 0) {
                    return this->GetValue();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_m->Differentiate(vars.at(0)->GetID());

                for (int i = 1; i < vars.size(); i++) {

                    temp = exp->Differentiate(vars.at(i)->GetID());
                    delete exp;
                    exp = temp;

                }
                T ret = exp->Evaluate();
                delete exp;

                return ret;
            }
        }

        /*!
         * Error of the derivative with respect to var0....var1 in order.
         * @param var0
         * @return numerical derivative
         */
        const T WRT_Error(const ADNumber<T> &var0) {
            //  Lock l(this->mutex_m);
            // Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            T ret = exp->PropagatedError();
            delete exp;
            return ret;
        }

        /*!
         * Error of the derivative with respect to var0....var2 in order.
         * @param var0
         * @param var1
         * @return numerical derivative
         */
        const T WRT_Error(const ADNumber<T> &var0, const ADNumber<T> &var1) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());

            T ret = temp->PropagatedError();
            delete temp;
            delete exp;


            return ret;
        }

        /*!
         * Error of the derivative with respect to var0....var3 in order.
         * @param var0
         * @param var1
         * @param var2
         * @return numerical derivative
         */
        const T WRT_Error(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            T ret = exp->PropagatedError();
            delete exp;

            return ret;
        }

        /*!
         * Error of the derivative with respect to var0....var4 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @return numerical derivative
         */
        const T WRT_Error(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            T ret = temp->Evaluate();
            delete temp;

            return ret;
        }

        /*!
         * Error of the derivative with respect to var0....var5 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @return numerical derivative
         */
        const T WRT_Error(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            T ret = exp->PropagatedError();
            delete exp;

            return ret;
        }

        /*!
         * Error of the derivative with respect to var0....var6 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @param var5
         * @return numerical derivative
         */
        const T WRT_Error(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4,
                const ADNumber<T> &var5) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            temp = exp->Differentiate(var5.GetID());
            delete exp;
            T ret = temp->Evaluate();
            delete temp;

            return ret;
        }

        /*!
         * Error of the derivative with respect to var0....var7 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @param var5
         * @param var6
         * @return numerical derivative
         */
        const T WRT_Error(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4,
                const ADNumber<T> &var5, const ADNumber<T> &var6) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            temp = exp->Differentiate(var5.GetID());
            delete exp;
            exp = temp->Differentiate(var6.GetID());
            delete temp;
            T ret = exp->PropagatedError();
            delete exp;

            return ret;
        }

        /*!
         * Error of the derivative with respect to var0....var8 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @param var5
         * @param var6
         * @param var7
         * @return numerical derivative
         */
        const T WRT_Error(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4,
                const ADNumber<T> &var5, const ADNumber<T> &var6, const ADNumber<T> &var7) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            temp = exp->Differentiate(var5.GetID());
            delete exp;
            exp = temp->Differentiate(var6.GetID());
            delete temp;
            temp = exp->Differentiate(var7.GetID());
            delete exp;
            T ret = temp->Evaluate();
            delete temp;

            return ret;
        }

        /*!
         * Error of the derivative with respect to var0....var9 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @param var5
         * @param var6
         * @param var7
         * @param var8
         * @return numerical derivative
         */
        const T WRT_Error(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4,
                const ADNumber<T> &var5, const ADNumber<T> &var6, const ADNumber<T> &var7,
                const ADNumber<T> &var8) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            temp = exp->Differentiate(var5.GetID());
            delete exp;
            exp = temp->Differentiate(var6.GetID());
            delete temp;
            temp = exp->Differentiate(var7.GetID());
            delete exp;
            exp = temp->Differentiate(var8.GetID());
            delete temp;
            T ret = exp->PropagatedError();
            delete exp;

            return ret;
        }

        /*!
         * Error of the derivative with respect to var0....var10 in order.
         * @param var0
         * @param var1
         * @param var2
         * @param var3
         * @param var4
         * @param var5
         * @param var6
         * @param var7
         * @param var8
         * @param var9
         * @return numerical derivative
         */
        const T WRT_Error(const ADNumber<T> &var0, const ADNumber<T> &var1,
                const ADNumber<T> &var2, const ADNumber<T> &var3, const ADNumber<T> &var4,
                const ADNumber<T> &var5, const ADNumber<T> &var6, const ADNumber<T> &var7,
                const ADNumber<T> &var8, const ADNumber<T> &var9) {
            //  Lock l(this->mutex_m);
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_m->Differentiate(var0.GetID());

            temp = exp->Differentiate(var1.GetID());
            delete exp;
            exp = temp->Differentiate(var2.GetID());
            delete temp;
            temp = exp->Differentiate(var3.GetID());
            delete exp;
            exp = temp->Differentiate(var4.GetID());
            delete temp;
            temp = exp->Differentiate(var5.GetID());
            delete exp;
            exp = temp->Differentiate(var6.GetID());
            delete temp;
            temp = exp->Differentiate(var7.GetID());
            delete exp;
            exp = temp->Differentiate(var8.GetID());
            delete temp;
            temp = exp->Differentiate(var9.GetID());
            delete exp;
            T ret = temp->Evaluate();
            delete temp;

            return ret;
        }

        /*!
         * Error of the derivative with respect to a vector of ADNumbers 
         * {var0,var1...varn} in order.
         * @param vars
         * @return numerical derivative
         */
        const T WRT_Error(const std::vector<ADNumber<T>*> &vars) {
            //  Lock l(this->mutex_m);
            if (this->expression_m == NULL) {
                return T(0);
            } else {
                if (vars.size() == 0) {
                    return this->GetValue();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_m->Differentiate(vars.at(0)->GetID());

                for (int i = 1; i < vars.size(); i++) {

                    temp = exp->Differentiate(vars.at(i)->GetID());
                    delete exp;
                    exp = temp;

                }
                T ret = exp->PropagatedError();
                delete exp;

                return ret;
            }
        }

        /*!
         * Return the nth order derivative.
         */
        const T Nth(const unsigned int &order) {
            //  Lock l(this->mutex_m);
            if (this->expression_m == NULL) {
                return T(0);
            } else {
                if (order == 0) {
                    return this->GetValue();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_m->Differentiate();

                for (size_t i = 1; i < order; i++) {

                    temp = exp->Differentiate();
                    delete exp;
                    exp = temp;

                }
                T ret = exp->Evaluate();
                delete exp;

                return ret;
            }

        }

        /*!
         * Return the nth order error of the derivative.
         */
        const T NthError(const unsigned int &order) {
            //  Lock l(this->mutex_m);
            if (this->expression_m == NULL) {
                return T(0);
            } else {
                if (order == 0) {
                    return this->expression_m->PropagatedError();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_m->Differentiate();

                for (size_t i = 1; i < order; i++) {

                    temp = exp->Differentiate();
                    delete exp;
                    exp = temp;

                }
                T ret = exp->PropagatedError();
                delete exp;

                return ret;
            }

        }

        /*!
         * Return the nth order partial derivative.
         */
        const T NthPartial(const ADNumber<T> &wrt, const unsigned int &order) {
            //  Lock l(this->mutex_m);
            if (!this->expression_m) {
                return T(0);
            } else {
                if (order == 0) {
                    return this->GetValue();
                }

                if (order == 1) {
                    return this->expression_m->EvaluateDerivative(wrt.GetID());
                }


                Expression<T> *temp;
                Expression<T> *exp = this->expression_m->Differentiate(wrt.GetID());


                if (order == 1) {
                    T ret = exp->Evaluate();
                    delete exp;

                    return ret;
                }


                size_t i;
                for (i = 1; i < order; i++) {

                    temp = exp->Differentiate(wrt.GetID());
                    delete exp;
                    exp = temp;

                }
                T ret = exp->Evaluate();
                delete exp;

                return ret;
            }

        }

        /*!
         * Return the nth order error of the partial derivative.
         */
        const T NthPartialError(const ADNumber<T> &wrt, unsigned int order) {
            //  Lock l(this->mutex_m);
            if (this->expression_m == NULL) {
                return T(0);
            } else {
                if (order == 0) {
                    return this->NthError(0);
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_m->Differentiate(wrt.GetID());

                for (int i = 1; i < order; i++) {

                    temp = exp->Differentiate(wrt.GetID());
                    delete exp;
                    exp = temp;

                }
                T ret = exp->PropagatedError();
                delete exp;

                return ret;
            }

        }

        /*!
         * Returns sqrt(sum a = A,B { (df/da)^2*u(a)^2}), where u(a) is 
         * round error.
         * 
         **/
        const T GetUncertainty() {
            //  Lock l(this->mutex_m);
            std::vector<uint32_t> vars;
            this->expression_m->VariableIds(vars);


            T temp = T(0);
            // T squared_epsilon = std::numeric_limits<T>::epsilon() * std::numeric_limits<T>::epsilon();
            T dif;
            for (size_t i = 0; i < vars.size(); i++) {
                Expression<T> *exp =
                        this->expression_m->Differentiate(vars.at(i));
                dif = exp->Evaluate();
                temp += dif * dif; //*squared_epsilon;
                delete exp;
            }

            return std::sqrt(temp);

        }



        //        
        //
        //        /*!
        //         * returns the gradient.
        //         * 
        //         * Harris, J and Stocker, H (1998). 
        //         * "Handbook of Mathematics and Computational Science", ''Springer'', 
        //         * page 513.
        //         * @return 
        //         */
        //        const T GetAbsoluteError() {
        //            return this->GetUncertainty();
        //        }
        //
        //        /*!
        //         * returns the gradient/value
        //         * Harris, J and Stocker, H (1998). 
        //         * "Handbook of Mathematics and Computational Science", ''Springer'', 
        //         * page 513.
        //         * @return 
        //         */
        //        const T GetRelativeError() {
        //            return GetAbsoluteError() / this->GetValue();
        //        }
        //
        //        /*!
        //         * returns the gradient/value * 100
        //         * Harris, J and Stocker, H (1998). 
        //         * "Handbook of Mathematics and Computational Science", ''Springer'', 
        //         * page 513.
        //         * @return 
        //         */
        //        const T GetPercentError() {
        //            return (GetAbsoluteError() / this->GetValue())*T(100);
        //        }
        //

        /*!
         * Return the nth order derivative as a std::string.
         */
        const std::string NthToString(const unsigned int &order) {
            //  Lock l(this->mutex_m);
            if (this->expression_m == NULL) {
                return "NA";
            } else {
                if (order == 0) {
                    return this->expression_m->ToString();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_m->Differentiate();

                for (int i = 1; i < order; i++) {
                    temp = exp->Differentiate();
                    delete exp;
                    exp = temp;

                }

                std::string ret = exp->ToString();
                delete exp;

                return ret;
            }

        }

        /*!
         * Return the nth order partial derivative as a std::string.
         */
        const std::string NthToCPPFunction(std::string name, const unsigned int &order) {
            //  Lock l(this->mutex_m);
            // std::cout << __func__ << ": Not yet Implemented!\n";
            //return "Not yet Implemented!\n";
            std::vector<std::string> args;

            std::string ret; // = this->expression_->ToString(args);


            if (this->expression_m == NULL) {
                return "NA";
            } else {
                if (order == 0) {
                    ret = this->expression_m->ToString(args);
                } else {

                    Expression<T> *temp;
                    Expression<T> *exp = this->expression_m->Differentiate();

                    for (int i = 1; i < order; i++) {
                        temp = exp->Differentiate();
                        delete exp;
                        exp = temp;

                    }

                    ret = exp->ToString(args);
                    delete exp;
                }
                std::stringstream ss;

                ss << "/*Machine generated by ADNumber::NthToCPPFunction*/\n";
                //                ss << "#ifdef AD_REAL\n";
                //                ss << "#undef AD_REAL\n";
                //                ss << "#endif\n\n";
                //                ss << "#define AD_REAL double\n\n\n";



                ss << "/*!\n"
                        " * Function " << name << ".\n */\n";

                ss << "template<class T> \nT " << name << "(";

                if (args.size() > 0) {
                    if (args.size() > 1) {

                        for (size_t i = 0; i < args.size() - 1; i++) {
                            ss << "T " << args.at(i) << ",";
                        }
                        ss << "T " << args.at(args.size() - 1) << "){\n";
                    } else {
                        ss << "T " << args.at(0) << "){\n";
                    }
                } else {
                    ss << "){";
                }

                ss << "\nT ret =" << ret << ";\n\nreturn ret;\n}";

                //   ss << "\n\n#undef AD_REAL\n";



                return ss.str();
            }

        }

        /*!
         * Derivative with respect to a vector of ADNumbers 
         * {var0,var1...varn} in order.
         * @param vars
         * @return numerical derivative
         */
        const std::string WRT_ToCPPFunction(std::string name, const std::vector<ADNumber<T> > &vars) {
            //  Lock l(this->mutex_m);
            // std::cout << __func__ << ": Not yet Implemented!\n";
            //return "Not yet Implemented!\n";
            std::vector<std::string> args;

            std::string ret; // = this->expression_->ToString(args);


            if (this->expression_m == NULL) {
                return "";
            } else {
                //                if (vars.size() == 0) {
                //                    ret = this->expression_m 
                //                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_m->Differentiate(vars.at(0).GetID());

                for (int i = 1; i < vars.size(); i++) {

                    temp = exp->Differentiate(vars.at(i).GetID());
                    delete exp;
                    exp = temp;

                }
                ret = exp->ToString();

                delete exp;

            }

            std::stringstream ss;

            ss << "/*Machine generated by ADNumber::NthToCPPFunction*/\n";
            //                ss << "#ifdef AD_REAL\n";
            //                ss << "#undef AD_REAL\n";
            //                ss << "#endif\n\n";
            //                ss << "#define AD_REAL double\n\n\n";



            ss << "/*!\n"
                    " * Function " << name << ".\n */\n";

            ss << "template<class T> \nT " << name << "(";

            if (args.size() > 0) {
                if (args.size() > 1) {

                    for (size_t i = 0; i < args.size() - 1; i++) {
                        ss << "T " << args.at(i) << ",";
                    }
                    ss << "T " << args.at(args.size() - 1) << "){\n";
                } else {
                    ss << "T " << args.at(0) << "){\n";
                }
            } else {
                ss << "){";
            }

            ss << "\nT ret =" << ret << ";\n\nreturn ret;\n}";

            //   ss << "\n\n#undef AD_REAL\n";



            return ss.str();



        }

        /*!
         * Return the nth order partial derivative as a std::string.
         */
        const std::string NthPartialToCPPFunction(std::string name, const ADNumber<T> &wrt, const unsigned int &order) {
            //  Lock l(this->mutex_m);
            // std::cout << __func__ << ": Not yet Implemented!\n";
            //return "Not yet Implemented!\n";
            std::vector<std::string> args;

            std::string ret; // = this->expression_->ToString(args);

            if (this->expression_m == NULL) {
                return "NA";
            } else {
                if (order == 0) {
                    ret = this->expression_m->ToString(args);
                } else {

                    Expression<T> *temp;
                    Expression<T> *exp = this->expression_m->Differentiate(wrt.GetID());

                    for (int i = 1; i < order; i++) {
                        temp = exp->Differentiate(wrt.GetID());
                        delete exp;
                        exp = temp;

                    }



                    ret = exp->ToString(args);
                    delete exp;

                }
                std::stringstream ss;

                ss << "/*Machine generated by ADNumber:NthPartialToCFuncton*/\n";




                ss << "/*!\n"
                        " * Function " << name << ".\n */\n";

                ss << "template<class T> \nT " << name << "(";

                if (args.size() > 0) {
                    if (args.size() > 1) {

                        for (size_t i = 0; i < args.size() - 1; i++) {
                            ss << "T " << args.at(i) << ",";
                        }
                        ss << "T " << args.at(args.size() - 1) << "){\n";
                    } else {
                        ss << "T " << args.at(0) << "){\n";
                    }
                } else {
                    ss << "){";
                }

                ss << "\nT ret =" << ret << ";\n\nreturn ret;\n}";





                return ss.str();
            }

        }

        /*!
         * Return the nth order partial derivative as a std::string.
         */
        const std::string NthPartialToString(const ADNumber<T> &wrt, const unsigned int &order) {
            //  Lock l(this->mutex_m);
            // std::cout << __func__ << ": Not yet Implemented!\n";
            //return "Not yet Implemented!\n";

            if (this->expression_m == NULL) {
                return "NA";
            } else {
                if (order == 0) {
                    return this->expression_m->ToString();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_m->Differentiate(wrt.GetID());

                for (int i = 1; i < order; i++) {
                    temp = exp->Differentiate(wrt.GetID());
                    delete exp;
                    exp = temp;

                }

                std::string ret = exp->ToString();
                delete exp;

                return ret;
            }

        }

        /*!
         * 
         * @return 
         */
        const std::string GetName() const {
            //  Lock l(this->mutex_m);
            return this->variableName_;
        }

        void SetName(const std::string &name) {
            //  Lock l(this->mutex_m);
            this->variableName_ = name;
        }

        /*
         * Return the unique identifier for this ADNumber.
         */
        const uint32_t GetID() const {
            //  Lock l(this->mutex_m);
            return this->id_m;
        }


        //expression tree, used for reverse mode calculation.
        Expression<T> *expression_m;

    protected:
        //computed value.
        T value_m;

        //forward computed derivative.(direct)
        T fderivative_m;


    private:
        //initialize this expression.

        void Initialize() {

            this->expression_m->SetValue(this->value_m);
            this->expression_m->SetId(this->id_m);
            this->expression_m->SetOp(VARIABLE);


        }


        //reverse computed derivative.holds the expression tree calculation to
        //avoid repeated evaluations of the expression tree.(adjoint)
        //T rderivative_;

        //reverse flag
        bool reverse_computed;

        //Default is x
        std::string variableName_;

        //unique id
        uint32_t id_m;

 


    public:
        //Friends
        // relational operators
        template<class TT> friend const int operator==(const ADNumber<TT>& lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const int operator!=(const ADNumber<TT>& lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const int operator<(const ADNumber<TT>& lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const int operator>(const ADNumber<TT>& lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const int operator<=(const ADNumber<TT>& lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const int operator>=(const ADNumber<TT>& lhs, const ADNumber<TT>& rhs);

        template<class TT> friend const int operator==(TT lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const int operator!=(TT lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const int operator<(TT lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const int operator>(TT lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const int operator<=(TT lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const int operator>=(TT lhs, const ADNumber<TT>& rhs);

        template<class TT> friend const int operator==(const ADNumber<TT>& lhs, TT rhs);
        template<class TT> friend const int operator!=(const ADNumber<TT>& lhs, TT rhs);
        template<class TT> friend const int operator<(const ADNumber<TT>& lhs, TT rhs);
        template<class TT> friend const int operator>(const ADNumber<TT>& lhs, TT rhs);
        template<class TT> friend const int operator<=(const ADNumber<TT>& lhs, TT rhs);
        template<class TT> friend const int operator>=(const ADNumber<TT>& lhs, TT rhs);


        // binary
        template<class TT> friend const ADNumber<TT> operator-(const ADNumber<TT>& lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const ADNumber<TT> operator/(const ADNumber<TT>& lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const ADNumber<TT> operator+(const ADNumber<TT>& lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const ADNumber<TT> operator*(const ADNumber<TT>& lhs, const ADNumber<TT>& rhs);


        template<class TT> friend const ADNumber<TT> operator-(TT lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const ADNumber<TT> operator/(TT lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const ADNumber<TT> operator+(TT lhs, const ADNumber<TT>& rhs);
        template<class TT> friend const ADNumber<TT> operator*(TT lhs, const ADNumber<TT>& rhs);


        template<class TT> friend const ADNumber<TT> operator-(const ADNumber<TT>& lhs, TT rhs);
        template<class TT> friend const ADNumber<TT> operator/(const ADNumber<TT>& lhs, TT rhs);
        template<class TT> friend const ADNumber<TT> operator+(const ADNumber<TT>& lhs, TT rhs);
        template<class TT> friend const ADNumber<TT> operator*(const ADNumber<TT>& lhs, TT rhs);


    };

    template<class T>
    ADNumber<T> Integrate(const ADNumber<T> &f, const ADNumber<T> &wrt) {
        ADNumber<T> ret;
        Expression<T>* tmp = f.expression_m->Integral(wrt.GetID());
        Expression<T>* forward = tmp->Differentiate();
        ret.Reset(f.expression_m->Evaluate(), forward->Evaluate());
        ret.expression_m = tmp;
        delete forward;
        return ret;
    }

    /*!
     * Equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator==(const ADNumber<T>& lhs, const ADNumber<T>& rhs) {

        return (lhs.GetValue() == rhs.GetValue());
    }

    /*!
     * Not equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator!=(const ADNumber<T>& lhs, const ADNumber<T>& rhs) {

        return (lhs.GetValue() != rhs.GetValue());
    }

    /*!
     * Less than comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator<(const ADNumber<T>& lhs, const ADNumber<T>& rhs) {

        return (lhs.GetValue() < rhs.GetValue());
    }

    /*!
     * Greater than comparison operator.
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator>(const ADNumber<T>& lhs, const ADNumber<T>& rhs) {

        return (lhs.GetValue() > rhs.GetValue());
    }

    /*!
     * Less than equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator<=(const ADNumber<T>& lhs, const ADNumber<T>& rhs) {

        return (lhs.GetValue() <= rhs.GetValue());
    }

    /*!
     * Greater than equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator>=(const ADNumber<T>& lhs, const ADNumber<T>& rhs) {

        return (lhs.GetValue() >= rhs.GetValue());
    }

    /*!
     * Equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator==(T lhs, const ADNumber<T>& rhs) {

        return (lhs == rhs.GetValue());
    }

    /*!
     * Not equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator!=(T lhs, const ADNumber<T>& rhs) {

        return (lhs != rhs.GetValue());
    }

    /*!
     * Less than comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator<(T lhs, const ADNumber<T>& rhs) {

        return (lhs < rhs.GetValue());
    }

    /*!
     * Greater than comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator>(T lhs, const ADNumber<T>& rhs) {

        return (lhs > rhs.GetValue());
    }

    /*!
     * Less than equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator<=(T lhs, const ADNumber<T>& rhs) {

        return (lhs <= rhs.GetValue());
    }

    /*!
     * Greater than equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator>=(T lhs, const ADNumber<T>& rhs) {

        return (lhs >= rhs.GetValue());
    }

    /*!
     * Equal to comparison operator.
     *
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator==(const ADNumber<T>& lhs, T rhs) {

        return (lhs.GetValue() == rhs);
    }

    /*!
     * Not equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator!=(const ADNumber<T>& lhs, T rhs) {

        return (lhs.GetValue() != rhs);
    }

    /*!
     * Less than comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator<(const ADNumber<T>& lhs, T rhs) {

        return (lhs.GetValue() < rhs);
    }

    /*!
     * Greater than comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator>(const ADNumber<T>& lhs, T rhs) {

        return (lhs.GetValue() > rhs);
    }

    /*!
     * Less than equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator<=(const ADNumber<T>& lhs, T rhs) {

        return (lhs.GetValue() <= rhs);
    }

    /*!
     * Greater than equal to comparison operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> inline const int operator>=(const ADNumber<T>& lhs, T rhs) {

        return (lhs.GetValue() >= rhs);
    }


    // binary

    /*!
     * Outside class subtraction operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator-(const ADNumber<T>& lhs, const ADNumber<T>& rhs) {
        ADNumber<T > ret(T(lhs.GetValue() - rhs.GetValue()),
                T(lhs.Forward() - rhs.Forward()));

        ret.expression_m->SetOp(MINUS);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Outside class addition operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator+(const ADNumber<T>& lhs, const ADNumber<T>& rhs) {
        ADNumber<T> ret(T(lhs.GetValue() + rhs.GetValue()),
                T(lhs.Forward() + rhs.Forward()));

        ret.expression_m->SetOp(PLUS);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Outside class division operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator/(const ADNumber<T>& lhs, const ADNumber<T>& rhs) {
        ADNumber<T > ret(T(lhs.GetValue() / rhs.GetValue()),
                T((rhs.GetValue() * lhs.Forward() - lhs.GetValue() * rhs.Forward())
                / (rhs.GetValue() * rhs.GetValue())));

        ret.expression_m->SetOp(DIVIDE);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Outside class multiplication operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator*(const ADNumber<T>& lhs, const ADNumber<T>& rhs) {
        ADNumber<T > ret(lhs.GetName(), T(lhs.GetValue() * rhs.GetValue()),
                T(lhs.GetValue() * rhs.Forward() + rhs.GetValue() * lhs.Forward()));

        ret.expression_m->SetOp(MULTIPLY);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Outside class subtraction operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator-(T lhs, const ADNumber<T>& rhs) {
        ADNumber<T > ret(T(lhs - rhs.GetValue()),
                T(0) - T(rhs.Forward()));

        Expression<T> *exp = new Expression<T > ();
        exp->SetValue(lhs);
        exp->SetOp(CONSTANT);

        ret.expression_m->SetOp(MINUS);
        ret.expression_m->SetLeft(exp);
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Outside class addition operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator+(T lhs, const ADNumber<T>& rhs) {
        ADNumber<T > ret(T(lhs + rhs.GetValue()),
                T(rhs.Forward()));

        Expression<T> *exp = new Expression<T > ();
        exp->SetValue(lhs);
        exp->SetOp(CONSTANT);

        ret.expression_m->SetOp(PLUS);
        ret.expression_m->SetLeft(exp);
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Outside class division operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator/(T lhs, const ADNumber<T>& rhs) {
        ADNumber<T > ret(T(lhs / rhs.GetValue()),
                T((rhs.GetValue() * 0 - lhs * rhs.Forward()) /
                (rhs.GetValue() * rhs.GetValue())));

        Expression<T> *exp = new Expression<T > ();
        exp->SetValue(lhs);
        exp->SetOp(CONSTANT);

        ret.expression_m->SetOp(DIVIDE);
        ret.expression_m->SetLeft(exp);
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;

    }

    /*!
     * Outside class multiplication operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator*(T lhs, const ADNumber<T>& rhs) {
        ADNumber<T > ret(T(lhs * rhs.GetValue()),
                T(lhs * rhs.Forward() + rhs.GetValue() * T(0)));

        Expression<T> *exp = new Expression<T > ();
        exp->SetValue(lhs);
        exp->SetOp(CONSTANT);

        ret.expression_m->SetOp(MULTIPLY);
        ret.expression_m->SetLeft(exp);
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Outside class subtraction operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator-(const ADNumber<T>& lhs, T rhs) {
        ADNumber<T > ret(T(lhs.GetValue() - rhs),
                T(lhs.Forward()));

        Expression<T> *exp = new Expression<T > ();
        exp->SetValue(rhs);
        exp->SetOp(CONSTANT);

        ret.expression_m->SetOp(MINUS);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(exp);

        return ret;
    }

    /*!
     * Outside class addition operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator+(const ADNumber<T>& lhs, T rhs) {
        ADNumber<T> ret(T(lhs.GetValue() + rhs),
                T(lhs.Forward()));

        Expression<T> *exp = new Expression<T > ();
        exp->SetValue(rhs);
        exp->SetOp(CONSTANT);

        ret.expression_m->SetOp(PLUS);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(exp);

        return ret;
    }

    /*!
     * Outside class division operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator/(const ADNumber<T>& lhs, T rhs) {
        ADNumber<T > ret(T(lhs.GetValue() / rhs),
                T((rhs * lhs.Forward() - lhs.GetValue() * 0) / (rhs * rhs)));

        Expression<T> *exp = new Expression<T > ();
        exp->SetValue(rhs);
        exp->SetOp(CONSTANT);

        ret.expression_m->SetOp(DIVIDE);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(exp);

        return ret;
    }

    /*!
     * Outside class multiplication operator.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> const ADNumber<T> operator*(const ADNumber<T>& lhs, T rhs) {
        ADNumber<T > ret(T(lhs.GetValue() * rhs),
                T(lhs.GetValue() * 0 + rhs * lhs.Forward()));

        Expression<T> *exp = new Expression<T > ();
        exp->SetValue(rhs);
        exp->SetOp(CONSTANT);

        ret.expression_m->SetOp(MULTIPLY);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(exp);

        return ret;
    }

    template<class T>
    std::ostream & operator<<(std::ostream &out, ADNumber<T> const &t) {
        out << t.GetValue();
        return out;
    }


} //ad
namespace std {

    //Math Overloads

    /*!
     * Returns the arc tangent of the ADNumber number val.
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> atan(const ad::ADNumber<T> &val) {
        ad::ADNumber<T> ret(atan(val.GetValue()),
                T(1.0) / (T(1.0) + pow(val.GetValue(), T(2))));

        ret.expression_m->SetOp(ad::ATAN);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    /*!
     * Compute arc tangent with two parameters. 
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> atan2(const ad::ADNumber<T> &lhs, const ad::ADNumber<T> &rhs) {

        T x = lhs.GetValue();
        T y = rhs.GetValue();
        T temp = x * x + y*y;
        ad::ADNumber<T> ret(atan2(x, y),
                (/*T(-1.0) * x*/y) / temp);

        ret.expression_m->SetOp(ad::ATAN2);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Compute arc tangent with two parameters.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> ad::ADNumber<T> atan2(T lhs, const ad::ADNumber<T> &rhs) {


        T x = lhs;
        T y = rhs.GetValue();
        T temp = x * x + y*y;

        ad::ADNumber<T> ret(atan2(x, y),
                (/*T(-1.0) * x*/y) / temp);

        ad::Expression<T> *exp = new ad::Expression<T > ();
        exp->SetValue(lhs);
        exp->SetOp(ad::CONSTANT);

        ret.expression_m->SetOp(ad::ATAN2);
        ret.expression_m->SetLeft(exp);
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Compute arc tangent with two parameters.
     *  
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> ad::ADNumber<T> atan2(const ad::ADNumber<T> &lhs, T rhs) {

        T x = lhs.GetValue();
        T y = rhs;
        T temp = x * x + y*y;
        ad::ADNumber<T> ret(atan2(x, y),
                y / temp);

        ad::Expression<T> *exp = new ad::Expression<T > ();
        exp->SetValue(rhs);
        exp->SetOp(ad::CONSTANT);

        ret.expression_m->SetOp(ad::ATAN2);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(exp);

        return ret;
    }

    /*!
     * Returns the cosine of the ADNumber number val.
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> cos(const ad::ADNumber<T> &val) {
        ad::ADNumber<T> ret(cos(val.GetValue()),
                T(-1) * sin(val.GetValue()));

        ret.expression_m->SetOp(ad::COS);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    /*!
     * Compute exponential function for val.
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> exp(const ad::ADNumber<T> &val) {
        ad::ADNumber<T> ret(exp(val.GetValue()),
                exp(val.GetValue()));

        ret.expression_m->SetOp(ad::EXP);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    template<class T> ad::ADNumber<T> mfexp(const ad::ADNumber<T> & x) {
        T b = T(60);
        if (x <= b && x >= T(-1) * b) {
            return std::exp(x);
        } else if (x > b) {
            return std::exp(b)*(T(1.) + T(2.) * (x - b)) / (T(1.) + x - b);
        } else {
            return std::exp(T(-1) * b)*(T(1.) - x - b) / (T(1.) + T(2.) * (T(-1) * x - b));
        }
    }

    /*!
     * Compute natural logarithm of val.
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> log(const ad::ADNumber<T> &val) {
        ad::ADNumber<T> ret(log(val.GetValue()), T(1.0) / val.GetValue());

        ret.expression_m->SetOp(ad::LOG);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    /*!
     * Compute natural common logarithm of val.
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> log10(const ad::ADNumber<T> &val) {
        ad::ADNumber<T> ret(log10(val.GetValue()),
                T(1.0) / (val.GetValue() * log(T(10.0))));



        ret.expression_m->SetOp(ad::LOG10);
        ret.expression_m->SetLeft(val.expression_m->Clone());


        return ret;
    }

    /*!
     * Raise to power.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> ad::ADNumber<T> pow(const ad::ADNumber<T> &lhs, const ad::ADNumber<T> &rhs) {
        ad::ADNumber<T> ret(pow(lhs.GetValue(), rhs.GetValue()),
                rhs.GetValue() * pow(lhs.GetValue(), rhs.GetValue() - T(1.0)));

        ret.expression_m->SetOp(ad::POW);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Raise to power.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> ad::ADNumber<T> pow(T lhs, const ad::ADNumber<T> & rhs) {
        ad::ADNumber<T> ret(pow(lhs, rhs.GetValue()),
                rhs.GetValue() * pow(lhs, rhs.GetValue() - T(1.0)));
        //   (lhs * (std::numeric_limits<T>::epsilon() / lhs)));

        ad::Expression<T> *exp = new ad::Expression<T > ();
        exp->SetValue(lhs);
        exp->SetOp(ad::CONSTANT);

        ret.expression_m->SetOp(ad::POW);
        ret.expression_m->SetLeft(exp);
        ret.expression_m->SetRight(rhs.expression_m->Clone());

        return ret;
    }

    /*!
     * Raise to power.
     * 
     * @param lhs
     * @param rhs
     * @return 
     */
    template<class T> ad::ADNumber<T> pow(const ad::ADNumber<T> &lhs, T rhs) {
        ad::ADNumber<T> ret(pow(lhs.GetValue(), rhs),
                rhs * pow(lhs.GetValue(), rhs - T(1.0)));

        ad::Expression<T> *exp = new ad::Expression<T > ();
        exp->SetValue(rhs);
        exp->SetOp(ad::CONSTANT);

        ret.expression_m->SetOp(ad::POW);
        ret.expression_m->SetLeft(lhs.expression_m->Clone());
        ret.expression_m->SetRight(exp);

        return ret;
    }

    /*!
     * Returns the sine of val.
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> sin(const ad::ADNumber<T> &val) {
        ad::ADNumber<T> ret(sin(val.GetValue()),
                cos(val.GetValue()));

        ret.expression_m->SetOp(ad::SIN);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    /*!
     * Compute square root val.
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> sqrt(const ad::ADNumber<T> &val) {
        T temp = sqrt(val.GetValue());
        ad::ADNumber<T> ret(temp,
                T(0.5) / temp);
        //        ad::Expression<T>* right = new ad::Expression<T > ();
        //        right->SetOp(ad::CONSTANT);
        //        right->SetValue(T(0.5));
        //
        //        //just use pow!!!
        //        ret.expression_m->SetOp(ad::POW);
        //        ret.expression_m->SetLeft(val.expression_m->Clone());
        //        ret.expression_m->SetRight(right);
        ret.expression_m->SetOp(ad::SQRT);
        ret.expression_m->SetLeft(val.expression_m->Clone());
        return ret;
    }

    /*!
     * Returns the tangent of val.
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> tan(const ad::ADNumber<T> &val) {
        T temp = cos(val.GetValue());
        ad::ADNumber<T> ret(tan(val.GetValue()),
                T(1.0) / (temp * temp));

        ret.expression_m->SetOp(ad::TAN);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    /*!
     * Returns the arc cosine of val.
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> acos(const ad::ADNumber<T> & val) {

        ad::ADNumber<T> ret(acos(val.GetValue()),
                T(-1.0) / sqrt(T(1.0) - pow(val.GetValue(), T(2))));

        ret.expression_m->SetOp(ad::ACOS);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    /*!
     * Returns the arc sine of val.
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> asin(const ad::ADNumber<T> &val) {
        ad::ADNumber<T> ret(asin(val.GetValue()), T(1.0) / sqrt(T(1.0) - pow(val.GetValue(), T(2))));

        ret.expression_m->SetOp(ad::ASIN);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    /*!
     * Returns the hyperbolic sin of val.
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> sinh(const ad::ADNumber<T> &val) {

        ad::ADNumber<T> ret(sinh(val.GetValue()),
                cosh(val.GetValue()));

        ret.expression_m->SetOp(ad::SINH);
        ret.expression_m->SetLeft(val.expression_m->Clone());
        return ret;
    }

    /*!
     * Returns the hyperbolic cosine of val.
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> cosh(const ad::ADNumber<T> &val) {
        ad::ADNumber<T> ret(cosh(val.GetValue()),
                sinh(val.GetValue()));

        ret.expression_m->SetOp(ad::COSH);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    /*!
     * Returns the hyperbolic tangent of val.
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> tanh(const ad::ADNumber<T> &val) {
        T temp = cosh(val.GetValue());
        ad::ADNumber<T> ret(std::tanh(val.GetValue()), (T(1) / temp)*(T(1) / temp));

        ret.expression_m->SetOp(ad::TANH);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    /*!
     * Compute absolute value of val.
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> fabs(const ad::ADNumber<T> &val) {

        ad::ADNumber<T> ret(fabs(val.GetValue()), fabs(val.Forward()));

        ret.expression_m->SetOp(ad::FABS);
        ret.expression_m->SetLeft(val.expression_m->Clone());

        return ret;
    }

    /*!
     * Round down value.
     * 
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> floor(const ad::ADNumber<T> &val) {
        ad::ADNumber<T> ret(floor(val.GetValue()), floor(val.Forward()));
        ret.expression_m->SetOp(ad::FLOOR);
        ret.expression_m->SetLeft(val.expression_m->Clone());
        return ret;
    }

    /*!
     * Solve an equation in the form:
     * lhs = rhs, for variable var;
     * @param lhs
     * @param rhs
     * @param var
     * @return 
     */
    template<class T> ad::ADNumber<T> solve(const ad::ADNumber<T> &lhs, const ad::ADNumber<T> &rhs, const ad::ADNumber<T> &var) {
        std::cout << "solve not yet implemented....\n";
        ad::ADNumber<T> ret;

        if (lhs.expression_m->HasID(var.GetID())) {
            std::cout << "left side contains var...\n";
        }

        if (rhs.expression_m->HasID(var.GetID())) {
            std::cout << "left side contains var...\n";
        }


        return ret;

    }


    //namespace std {

    template<class T>
    class numeric_limits<ad::ADNumber<T> > : public numeric_limits<T> {
    };





}

typedef ad::ADNumber<double> addouble;
typedef ad::ADNumber<double> adfloat;

typedef ad::ADNumber<double> dvar;
typedef ad::ADNumber<double> fvar;

#endif	/* ADNUMBER_HPP */
