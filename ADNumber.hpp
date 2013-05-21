#ifndef ADDNUMBER_HPP
#define	ADDNUMBER_HPP

/*!
 *  Software to compute derivatives. Support for forward, reverse, partial, 
 *  nth, and nth partial derivatives.
 */

/*!
 *   This library is dual-licensed: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License version 3 as 
 *   published by the Free Software Foundation. For the terms of this 
 *   license, see licenses/gpl_v3.txt or <http://www.gnu.org/licenses/>.
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
 * Email:  msupernaw@gmail.com
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

//#define USE_MEMORY_POOL

//#include <complex>
#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <limits>
#include <stdint.h>

#ifdef USE_MEMORY_POOL
#include "memory_pool.hpp"
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

        /*!
         * Constructor.
         * 
         * @param id
         * @param op
         * @param value
         * @param left
         * @param right
         */
        Expression(unsigned int id, Operation op,
                T value,
                Expression<T>* left,
                Expression<T>* right) :
        id_(id),
        op_(op),
        value_(value),
        left_(left),
        right_(right),
        epsilon_m(std::numeric_limits<T>::epsilon()) {

        }

        /*!
         * Default constructor.
         */
        Expression() :
        left_(NULL),
        right_(NULL),
        value_(T(1.0)),
        id_(0),
        op_(VARIABLE),
        epsilon_m(std::numeric_limits<T>::epsilon()) {

        }
        
        /*!
        *Copy Constructor.
        */
        Expression(const Expression &orig) :
        id_(orig.id_),
        op_(orig.op_),
        value_(orig.value_),
        left_(NULL),
        right_(NULL),
        epsilon_m(std::numeric_limits<T>::epsilon()) {

            if (orig.left_ != NULL) {
                this->left_ = orig.left_.Clone();
            }
            if (orig.right_ != NULL) {
                this->right_ = orig.right_.Clone();
            }

            // this = orig.Clone();
        }

        /*!
         * Destructor.
         */
        virtual ~Expression() {

            if (this->left_) {
                delete this->left_;
            }

            if (this->right_) {
                delete this->right_;
            }
        }
#ifdef USE_MEMORY_POOL

        void* operator new (size_t size) throw (std::bad_alloc) {
            assert(size == sizeof (Expression));
            void* ptr = MemoryPool < sizeof (Expression)>::instance().allocate();
            return ptr;
        }

        void operator delete (void* ptr)throw () {

            MemoryPool < sizeof (Expression)>::instance().deallocate((void*) ptr);
        }
#endif

        /*!
         * Create a clone of this expression. The same as using
         * copy constructor.
         * @return 
         */
        Expression<T>* Clone() {
            Expression<T> *exp = new Expression<T > (this->id_, this->op_, this->value_, NULL, NULL);
            exp->op_ = this->op_;
            exp->id_ = this->id_;
            exp->value_ = this->value_;
            //
            //            if (this->right_ && this->left_) {
            //                return new Expression<T > (this->id_, this->op_, this->value_, this->left_->Clone(), this->right_->Clone());
            //            }

            if (this->left_) {
                //  return new Expression<T > (this->id_, this->op_, this->value_, this->left_->Clone(), NULL);
                exp->left_ = this->left_->Clone();
            }
            if (this->right_) {
                // return new Expression<T > (this->id_, this->op_, this->value_, NULL, this->right_->Clone());
                exp->right_ = this->right_->Clone();
            }



            return exp; //new Expression<T > (this->id_, this->op_, this->value_, NULL, NULL); //exp;
        }

        /*!
         * Evaluate this expression. 
         * @return 
         */
        const T Evaluate() const {
            switch (op_) {
                case CONSTANT:
                    return this->value_;
                case VARIABLE:
                    return this->value_;
                case MINUS:
                    return (this->left_->Evaluate() - this->right_->Evaluate());
                case PLUS:
                    return (this->left_->Evaluate() + this->right_->Evaluate());
                case DIVIDE:
                    return (this->left_->Evaluate() / this->right_->Evaluate());
                case MULTIPLY:
                    return (this->left_->Evaluate() * this->right_->Evaluate());
                case SIN:
                    return sin(this->left_->Evaluate());
                case COS:
                    return cos(this->left_->Evaluate());
                case TAN:
                    return tan(this->left_->Evaluate());
                case ASIN:
                    return asin(this->left_->Evaluate());
                case ACOS:
                    return acos(this->left_->Evaluate());
                case ATAN:
                    return atan(this->left_->Evaluate());
                case ATAN2:
                    return atan2(this->left_->Evaluate(), this->right_->Evaluate());
                    //                case ATAN3:
                    //                    break;
                    //                case ATAN4:
                    //                    break;
                case SQRT:
                    return sqrt(this->left_->Evaluate());
                case POW:
                    return pow(this->left_->Evaluate(), this->right_->Evaluate());
                    //                case POW1:
                    //                    break;
                    //                case POW2:
                    //                    break;
                case LOG:
                    return log(this->left_->Evaluate());
                case LOG10:
                    return log10(this->left_->Evaluate());
                case EXP:
                    return exp(this->left_->Evaluate());
                case SINH:
                    return sinh(this->left_->Evaluate());
                case COSH:
                    return cosh(this->left_->Evaluate());
                case TANH:
                    return tanh(this->left_->Evaluate());
                case FABS:
                    return fabs(this->left_->Evaluate());
                case ABS:
                    return abs(this->left_->Evaluate());
                case FLOOR:
                    return floor(this->left_->Evaluate());
                case NONE:
                    return this->value_;
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

            switch (op_) {


                case CONSTANT:
                    return this->epsilon_m;
                case VARIABLE:
                    return this->epsilon_m;
                case MINUS:
                    return (this->left_->PropagatedError() * this->left_->PropagatedError() + this->right_->PropagatedError() * this->right_->PropagatedError());
                case PLUS:
                    return (this->left_->PropagatedError() * this->left_->PropagatedError() + this->right_->PropagatedError() * this->right_->PropagatedError());
                case DIVIDE:

                    if (this->left_) {
                        a = this->left_->Evaluate();
                    } else {
                        a = T(0);
                    }
                    err_a = this->left_->PropagatedError();
                    if (this->right_) {
                        b = this->right_->Evaluate();
                    } else {
                        b = T(0);
                    }
                    err_b = this->right_->PropagatedError();
                    return std::sqrt((err_a * err_a) / (b * b) + (a * a) *(err_b * err_b)*(b * b * b * b));

                case MULTIPLY:

                    if (this->left_) {
                        a = this->left_->Evaluate();
                    } else {
                        a = T(0);
                    }
                    err_a = this->left_->PropagatedError();
                    if (this->right_) {
                        b = this->right_->Evaluate();
                    } else {
                        b = T(0);
                    }
                    err_b = this->right_->PropagatedError();

                    return std::sqrt((b * b)* (err_a * err_a)+(a * a)*(err_b * err_b));

                case SIN:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();

                    return std::fabs(std::cos(a) * err_a);
                case COS:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();

                    return std::fabs(std::sin(a) * err_a);
                case TAN:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();
                    return std::fabs(err_a / std::pow(std::cos(a), T(2.0)));
                case ASIN:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();
                    return (err_a / sqrt(T(1.0) - pow(a, T(2.0))));
                case ACOS:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();
                    return (err_a / sqrt(T(1.0) + pow(a, 2.0)));
                case ATAN:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();
                    return (err_a / sqrt(T(1.0) + pow(a, 2.0)));
                case ATAN2:

                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();
                    b = this->right_->Evaluate();
                    err_b = this->right_->PropagatedError();

                    temp_ = fabs(a / b) *
                            sqrt(pow(err_a / a, T(2.0)) +
                            pow(err_b / b, T(2.0)));
                    error = fabs(temp_ / (T(1.0) + pow(a / b, T(2.0))));

                    return error;
                case ATAN3:
                    break;
                case ATAN4:
                case SQRT:

                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();

                    return std::fabs(err_a / T(2) * sqrt(a));
                case POW:

                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();
                    b = this->right_->Evaluate();
                    err_b = this->right_->PropagatedError();

                    return sqrt((b * b) * std::pow(a, (b - T(1)))*(err_a * err_a) + std::log(a) * std::log(a) * std::pow(a, b)*(err_b * err_b));
                case POW1:
                    break;
                case POW2:
                    break;
                case LOG:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();
                    return std::fabs(err_a / a);
                case LOG10:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();
                    return std::fabs(err_a / a * std::log10(T(10.0)));
                case EXP:
                    if (this->left_) {
                        a = this->left_->Evaluate();
                        b = this->left_->Evaluate();
                    } else {
                        a = std::numeric_limits<T>::epsilon();
                        b = std::numeric_limits<T>::epsilon();
                    }

                    return std::fabs(std::exp(a) * b);
                case SINH:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();

                    return std::fabs(std::cosh(a) * err_a);
                case COSH:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();

                    return std::fabs(std::sinh(a) * err_a);
                case TANH:
                    a = this->left_->Evaluate();
                    err_a = this->left_->PropagatedError();
                    return std::fabs(err_a / std::pow(std::cosh(a), T(2.0)));
                case FABS:
                    return this->left_->PropagatedError();
                case ABS:
                    return this->left_->PropagatedError();
                case FLOOR:
                    return this->left_->PropagatedError();
                case NONE:
                    return this->left_->PropagatedError();
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

            Expression<T>* temp;
            switch (op_) {

                case CONSTANT:
                    //f(x) = C
                    //f'(x) = 0

                    ret->op_ = CONSTANT;
                    ret->value_ = T(0.0);

                    //                    //ret->Simplify();
                    return ret;

                case VARIABLE:
                    //f(x) = x
                    //f'(x) = 1

                    ret->op_ = CONSTANT;
                    ret->value_ = T(1.0);

                    //ret->Simplify();
                    return ret;
                case MINUS:
                    //f(x) = g(x) - h(x)
                    //f'(x) = g'(x) - h'(x)

                    ret->op_ = MINUS;
                    if (this->left_ != NULL) {
                        ret->left_ = this->left_->Differentiate();

                    }

                    if (this->right_ != NULL) {
                        ret->right_ = this->right_->Differentiate();
                    }
                    //ret->Simplify();
                    return ret;
                case PLUS:
                    //f(x) = g(x) + h(x)
                    //f'(x) = g'(x) + h'(x)

                    ret->op_ = PLUS;
                    if (this->left_ != NULL) {
                        ret->left_ = this->left_->Differentiate();
                    }

                    if (this->right_ != NULL) {
                        ret->right_ = this->right_->Differentiate();
                    }

                    //ret->Simplify();
                    return ret;
                case DIVIDE:
                    //f(x) = g(x)/h(x);
                    //f'(x) = (g'(x)h(x) - g(x)h'(x))/h(x)^2

                    ret->op_ = DIVIDE;

                    ret->left_ = new Expression<T > (); //g'(x)h(x) - g(x)h'(x)
                    ret->left_->op_ = MINUS;


                    ret->left_->left_ = new Expression<T > (); //g'(x)h(x)
                    ret->left_->left_->op_ = MULTIPLY;
                    if (this->left_ != NULL) {
                        ret->left_->left_->left_ = this->left_->Differentiate();
                    }
                    ret->left_->left_->right_ = this->right_->Clone();

                    ret->left_->right_ = new Expression<T > (); //g(x)h'(x)
                    ret->left_->right_->op_ = MULTIPLY;
                    ret->left_->right_->left_ = this->left_->Clone();
                    if (this->right_ != NULL) {
                        ret->left_->right_->right_ = this->right_->Differentiate();
                    }


                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = MULTIPLY;
                    ret->right_->left_ = this->right_->Clone();
                    ret->right_->right_ = this->right_->Clone();

                    //ret->Simplify();
                    return ret;

                case MULTIPLY:
                    //f(x) = g(x)h(x);
                    //f'(x) = g'(x)h(x) + g(x)h'(x)

                    if (this->left_->op_ == CONSTANT
                            && this->right_->op_ != CONSTANT) {
                        ret->op_ = MULTIPLY;

                        ret->left_ = this->left_->Clone();
                        ret->right_ = this->right_->Differentiate();


                    } else if (this->right_->op_ == CONSTANT
                            && this->left_->op_ != CONSTANT) {
                        ret->op_ = MULTIPLY;

                        ret->left_ = this->left_->Differentiate();
                        ret->right_ = this->right_->Clone();

                    } else {



                        ret->op_ = PLUS;

                        ret->left_ = new Expression<T > ();
                        ret->left_->op_ = MULTIPLY;

                        ret->left_->right_ = this->right_->Clone();

                        if (this->right_ != NULL) {
                            ret->left_->left_ = this->left_->Differentiate();
                        }
                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = MULTIPLY;

                        ret->right_->left_ = this->left_->Clone();
                        if (this->left_ != NULL) {
                            ret->right_->right_ = this->right_->Differentiate();
                        }



                    }
                    //ret->Simplify();
                    return ret;

                case SIN:
                    //f'(x) = cos(x)

                    ret->op_ = MULTIPLY;

                    ret->left_ = this->left_->Differentiate();

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = COS;
                    ret->right_->left_ = this->left_->Clone();

                    //ret->Simplify();
                    return ret;

                case COS:
                    //f'(x) = -sin(x)


                    ret->op_ = MULTIPLY;


                    ret->left_ = new Expression<T > ();
                    ret->left_->op_ = MULTIPLY;
                    ret->left_->left_ = new Expression<T>;
                    ret->left_->left_->op_ = CONSTANT;
                    ret->left_->left_->value_ = T(-1.0);
                    ret->left_->right_ = this->left_->Differentiate();


                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = SIN;
                    ret->right_->left_ = this->left_->Clone();

                    //ret->Simplify();
                    return ret;
                case TAN:
                    //f(x) = tan(x)
                    //f'(x) = (1/cos(x))(1/cos(x))

                    //easier to use the identity.
                    temp = new Expression<T > ();
                    temp->op_ = DIVIDE;

                    temp->left_ = new Expression<T>;
                    temp->left_->op_ = SIN;
                    temp->left_->left_ = this->left_->Clone();

                    temp->right_ = new Expression<T>;
                    temp->right_->op_ = COS;
                    temp->right_->left_ = this->left_->Clone();

                    ret = temp->Differentiate();



                    //                    
                    //                    ret->op_ = MULTIPLY;
                    //
                    //                    ret->left_ = new Expression<T > ();
                    //                    ret->left_->op_ = DIVIDE;
                    //
                    //                    ret->left_->left_ = new Expression<T > ();
                    //                    ret->left_->left_->op_ = CONSTANT;
                    //                    ret->left_->left_->value_ = T(1.0);
                    //
                    //                    ret->left_->right_ = new Expression<T > ();
                    //                    ret->left_->right_->op_ = COS;
                    //                    ret->left_->right_->left_ = this->left_->Clone();
                    //
                    //                    ret->right_ = ret->left_->Clone();

                    //ret->Simplify();
                    return ret;
                case ASIN:
                    //f(x) = asin(x)
                    //f'(x) = 1/(2 sqrt(1-x^2)= 1/(pow((1-pow(x,2)),0.5)

                    ret->op_ = DIVIDE;

                    ret->left_ = new Expression<T > ();
                    ret->left_->op_ = CONSTANT;
                    ret->left_ ->value_ = T(1.0);

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = POW;

                    ret->right_->left_ = new Expression<T > ();
                    ret->right_->left_->op_ = MINUS;

                    ret->right_->left_->left_ = new Expression<T > ();
                    ret->right_->left_->left_->op_ = CONSTANT;
                    ret->right_->left_->left_->value_ = T(1.0);

                    ret->right_->left_->right_ = new Expression<T > ();
                    ret->right_->left_->right_->op_ = POW;
                    ret->right_->left_->right_->left_ = this->left_->Clone();

                    ret->right_->left_->right_->right_ = new Expression<T > ();
                    ret->right_->left_->right_->right_->op_ = CONSTANT;
                    ret->right_->left_->right_->right_->value_ = T(2.0);

                    ret->right_->right_ = new Expression<T > ();
                    ret->right_->right_->op_ = CONSTANT;
                    ret->right_->right_->value_ = T(0.5);
                    //ret->Simplify();
                    return ret;
                case ACOS:
                    //f(x) = acos(x)
                    //f'(x) = -1/(sqrt(1-x^2) = -1/(pow((1-pow(x,2)),0.5)
                    //-1/sqrt(1-x^2)

                    ret->op_ = DIVIDE;

                    ret->left_ = new Expression<T > ();
                    ret->left_->op_ = CONSTANT;
                    ret->left_ ->value_ = T(-1.0);

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = POW;

                    ret->right_->left_ = new Expression<T > ();
                    ret->right_->left_->op_ = MINUS;

                    ret->right_->left_->left_ = new Expression<T > ();
                    ret->right_->left_->left_->op_ = CONSTANT;
                    ret->right_->left_->left_->value_ = T(1.0);

                    ret->right_->left_->right_ = new Expression<T > ();
                    ret->right_->left_->right_->op_ = POW;
                    ret->right_->left_->right_->left_ = this->left_->Clone();

                    ret->right_->left_->right_->right_ = new Expression<T > ();
                    ret->right_->left_->right_->right_->op_ = CONSTANT;
                    ret->right_->left_->right_->right_->value_ = T(2.0);

                    ret->right_->right_ = new Expression<T > ();
                    ret->right_->right_->op_ = CONSTANT;
                    ret->right_->right_->value_ = T(0.5);
                    //ret->Simplify();
                    return ret;
                case ATAN:
                    //f(x) = atan(x)
                    //f'(x) 1/(x^2+1)

                    ret->op_ = DIVIDE;
                    ret->left_ = new Expression<T > ();
                    ret->left_->op_ = MULTIPLY;
                    ret->left_ ->left_ = new Expression<T > ();
                    ret->left_->left_->op_ = CONSTANT;
                    ret->left_->left_->value_ = T(1.0);
                    ret->left_->right_ = this->left_->Differentiate();

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = PLUS;

                    ret->right_->left_ = new Expression<T > ();
                    ret->right_->left_->op_ = MULTIPLY;
                    ret->right_->left_->left_ = this->left_->Clone();
                    ret->right_->left_->right_ = this->left_->Clone();


                    ret->right_->right_ = new Expression<T > ();
                    ret->right_->right_->op_ = CONSTANT;
                    ret->right_->right_->value_ = T(1.0);

                    //ret->Simplify();
                    return ret;
                case ATAN2:
                    //f(x) = atan2(x,y)
                    //f'(x) y/(x^2+y^2)

                    ret->op_ = DIVIDE;
                    ret->left_ = new Expression<T > ();
                    ret->left_->op_ = MULTIPLY;
                    ret->left_ ->left_ = this->right_->Clone();
                    ret->left_->right_ = this->left_->Differentiate();


                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = PLUS;

                    ret->right_->left_ = new Expression<T > ();
                    ret->right_->left_->op_ = MULTIPLY;
                    ret->right_->left_->left_ = this->left_->Clone();
                    ret->right_->left_->right_ = this->left_->Clone();


                    ret->right_->right_ = new Expression<T > ();
                    ret->right_->right_->op_ = MULTIPLY;
                    ret->right_->right_->left_ = this->right_->Clone();
                    ret->right_->right_->right_ = this->right_->Clone();

                    //ret->Simplify();
                    return ret;

                    //  case ATAN4:
                case SQRT:
                    //f(x) = sqrt(x)
                    //f'(x) = .5/sqrt(x)

                    ret->op_ = DIVIDE;
                    ret->left_ = new Expression<T > ();
                    ret->left_->op_ = CONSTANT;
                    ret->left_->value_ = T(0.5);

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = SQRT;
                    ret->right_->left_ = this->left_->Clone();


                    //ret->Simplify();
                    return ret;
                case POW:
                    //f(x) =  x^y
                    //f'(x) = yx^y-1

                    ret->op_ = MULTIPLY;

                    ret->left_ = new Expression<T > ();
                    ret->left_->op_ = MULTIPLY;
                    ret->left_->left_ = this->left_->Differentiate();
                    ret->left_->right_ = this->right_->Clone();



                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = POW;


                    ret->right_->left_ = this->left_->Clone();


                    ret->right_->right_ = new Expression<T > ();
                    ret->right_->right_->op_ = MINUS;
                    ret->right_->right_->left_ = this->right_->Clone();

                    ret->right_->right_->right_ = new Expression<T > ();
                    ret->right_->right_->right_->op_ = CONSTANT;
                    ret->right_->right_->right_->value_ = T(1.0);

                    //ret->Simplify();

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

                    ret->op_ = DIVIDE;
                    ret->left_ = new Expression<T > ();
                    ret->left_->op_ = CONSTANT;
                    ret->left_->value_ = T(1.0);

                    ret->right_ = this->left_->Clone();


                    //ret->Simplify();
                    return ret;
                case LOG10:
                    //f(x) = log10(x)
                    //f'(x) = 1/(xlog(10))

                    ret->op_ = DIVIDE;
                    ret->left_ = new Expression<T > ();
                    ret->left_->op_ = CONSTANT;
                    ret->left_->value_ = T(1.0);

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = MULTIPLY;

                    ret->right_->left_ = this->left_->Clone();

                    ret->right_->right_ = new Expression<T > ();
                    ret->right_->right_->op_ = CONSTANT;
                    ret->right_->right_->value_ = log(T(10.0));
                    /*
                    ret->right_->right_->left_ = new Expression<T > ();
                    ret->right_->right_->left_->op_ = CONSTANT;
                    ret->right_->right_->left_->value_ = T(10.0);
                     */
                    //ret->Simplify();
                    return ret;
                case EXP:
                    //f(x) = e^x
                    //f'(x) =e^x

                    ret->op_ = MULTIPLY;
                    ret->left_ = this->left_->Differentiate();

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = EXP;
                    ret->right_->left_ = this->left_->Clone();

                    //ret->Simplify();
                    return ret;
                case SINH:
                    //f(x) = sinh(x)
                    //f'(x) = cosh(x)
                    ret->op_ = MULTIPLY;

                    ret->left_ = this->left_->Differentiate();

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = COSH;
                    ret->right_->left_ = this->left_->Clone();

                    //ret->Simplify();
                    return ret;
                case COSH:
                       ret->op_ = MULTIPLY;

                    ret->left_ = this->left_->Differentiate();

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = SINH;
                    ret->right_->left_ = this->left_->Clone();


                    //ret->Simplify();
                    return ret;
                case TANH:
                    //f(x) = tanh(x)
                    //f'(x) =1- tanh(x)*tanh(x)

                    ret->op_ = MINUS;

                    ret->left_ = new Expression<T > ();
                    ret->left_->op_ = CONSTANT;
                    ret->left_->value_ = T(1.0);


                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = MULTIPLY;

                    ret->right_->left_ = this->Clone();
                    ret->right_->right_ = this->Clone();

                    //ret->Simplify();
                    return ret;

                case FABS:

                    ret->op_ = DIVIDE;
                    ret->left_ = this->left_->Clone();

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = FABS;
                    ret->right_->left_ = this->left_->Clone();

                    //ret->Simplify();

                    return ret;
                case ABS:
#warning needs review
                    ret->op_ = DIVIDE;
                    ret->left_ = this->left_->Clone();

                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = ABS;
                    ret->right_->left_ = this->left_->Clone();

                    //ret->Simplify();
                    return ret;

                case FLOOR:

                    ret->op_ = FLOOR;
                    ret->left_ = this->left_->Clone();

                    //ret->Simplify();
                    return ret;
                case NONE://shouldn't happen.
                    return this->Clone();

                default:
                    return NULL;
            }
            return NULL;
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


            switch (op_) {

                case CONSTANT:
                    //f(x) = C
                    //f'(x) = 0

                    ret->op_ = CONSTANT;
                    ret->value_ = this->value_;

                    //ret->Simplify();
                    return ret;

                case VARIABLE:
                    if (this->id_ == id) {
                        //f(x) = x
                        //f'(x) = 1

                        ret->op_ = CONSTANT;
                        ret->value_ = T(1.0);

                        //ret->Simplify();
                        return ret;
                    } else {//constant
                        //f(x) = C
                        //f'(x) = 0
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);
                        return ret;
                    }
                case MINUS:

                    //f(x) = g(x) - h(x)
                    //f'(x) = g'(x) - h'(x)

                    ret->op_ = MINUS;
                    if (this->left_) {
                        ret->left_ = this->left_->Differentiate(id);

                    }

                    if (this->right_) {
                        ret->right_ = this->right_->Differentiate(id);
                    }
                    //ret->Simplify();
                    return ret;

                case PLUS:

                    //f(x) = g(x) + h(x)
                    //f'(x) = g'(x) + h'(x)

                    ret->op_ = PLUS;
                    if (this->left_) {
                        ret->left_ = this->left_->Differentiate(id);
                    }

                    if (this->right_) {
                        ret->right_ = this->right_->Differentiate(id);
                    }

                    //ret->Simplify();
                    return ret;

                case DIVIDE:

                    //f(x) = g(x)/h(x);
                    //f'(x) = (g'(x)h(x) - g(x)h'(x))/h(x)^2


                    ret->op_ = DIVIDE;

                    ret->left_ = new Expression<T > (); //g'(x)h(x) - g(x)h'(x)
                    ret->left_->op_ = MINUS;


                    ret->left_->left_ = new Expression<T > (); //g'(x)h(x)
                    ret->left_->left_->op_ = MULTIPLY;
                    if (this->left_) {
                        ret->left_->left_->left_ = this->left_->Differentiate(id);
                    }
                    ret->left_->left_->right_ = this->right_->Clone();

                    ret->left_->right_ = new Expression<T > (); //g(x)h'(x)
                    ret->left_->right_->op_ = MULTIPLY;
                    ret->left_->right_->left_ = this->left_->Clone();
                    if (this->right_) {
                        ret->left_->right_->right_ = this->right_->Differentiate(id);
                    }


                    ret->right_ = new Expression<T > ();
                    ret->right_->op_ = MULTIPLY;
                    ret->right_->left_ = this->right_->Clone();
                    ret->right_->right_ = this->right_->Clone(); /*= new Expression<T > ();
                ret->right_->right_->op_ = CONSTANT;
                ret->right_->right_->value_ = T(2.0);

*/

                    //ret->Simplify();
                    return ret;

                case MULTIPLY:
                    //f(x) = g(x)h(x);
                    //f'(x) = g'(x)h(x) + g(x)h'(x)

                    if (this->left_->op_ == CONSTANT
                            && this->right_->op_ != CONSTANT) {
                        ret->op_ = MULTIPLY;

                        ret->left_ = this->left_->Clone();
                        ret->right_ = this->right_->Differentiate(id);


                    } else if (this->right_->op_ == CONSTANT
                            && this->left_->op_ != CONSTANT) {
                        ret->op_ = MULTIPLY;

                        ret->left_ = this->left_->Differentiate(id);
                        ret->right_ = this->right_->Clone();

                    } else {



                        ret->op_ = PLUS;

                        ret->left_ = new Expression<T > ();
                        ret->left_->op_ = MULTIPLY;

                        ret->left_->right_ = this->right_->Clone();

                        if (this->right_ != NULL) {
                            ret->left_->left_ = this->left_->Differentiate(id);
                        }
                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = MULTIPLY;

                        ret->right_->left_ = this->left_->Clone();
                        if (this->left_ != NULL) {
                            ret->right_->right_ = this->right_->Differentiate(id);
                        }


                        //ret->Simplify();
                    }
                    return ret;

                case SIN:

                    if (this->left_->HasID(id)) {
                        //f'(x) = cos(x)

                        ret->op_ = MULTIPLY;
                        ret->left_ = this->left_->Differentiate(id);
                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = COS;
                        ret->right_->left_ = this->left_->Clone();
                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }

                case COS:
                    if (this->left_->HasID(id)) {
                        //f'(x) = -sin(x)

                        ret->op_ = MULTIPLY;


                        ret->left_ = this->left_->Differentiate(id);
                        ret->right_ = new Expression<T > ();

                        ret->right_->op_ = MULTIPLY;
                        ret->right_->left_ = new Expression<T > ();
                        ret->right_->left_->op_ = CONSTANT;
                        ret->right_->left_->value_ = T(-1.0);

                        ret->right_->right_ = new Expression<T > ();
                        ret->right_->right_->op_ = SIN;
                        ret->right_->right_->left_ = this->left_->Clone();

                        //ret->Simplify();
                        return ret;

                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case TAN:
                    if (this->left_->HasID(id)) {
                        //f'(x) = 1/cos(x)

                        ret->op_ = MULTIPLY;
                        ret->left_ = this->left_->Differentiate(id);


                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = MULTIPLY;

                        ret->right_->left_ = new Expression<T > ();
                        ret->right_->left_->op_ = DIVIDE;


                        ret->right_->left_->left_ = new Expression<T > ();
                        ret->right_->left_->left_->op_ = CONSTANT;
                        ret->right_->left_->left_->value_ = T(1.0);


                        ret->right_->left_->right_ = new Expression<T > ();
                        ret->right_->left_->right_->op_ = COS;
                        ret->right_->left_->right_->left_ = this->left_->Clone();


                        ret->right_->right_ = new Expression<T > ();
                        ret->right_->right_->op_ = DIVIDE;


                        ret->right_->right_->left_ = new Expression<T > ();
                        ret->right_->right_->left_->op_ = CONSTANT;
                        ret->right_->right_->left_->value_ = T(1.0);


                        ret->right_->right_->right_ = new Expression<T > ();
                        ret->right_->right_->right_->op_ = COS;
                        ret->right_->right_->right_->left_ = this->left_->Clone();




                        //
                        //                ret->right_ = new Expression<T > ();
                        //                ret->right_->op_ = DIVIDE;
                        //
                        //                ret->right_->left_ = new Expression<T > ();
                        //                ret->right_->left_->op_ = CONSTANT;
                        //                ret->right_->left_->value_ = T(1.0);
                        //
                        //                ret->right_->right_ = new Expression<T > ();
                        //                ret->right_->right_->op_ = COS;
                        //                ret->right_->right_->left_ = this->left_->Clone();



                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case ASIN:

                    if (this->left_->HasID(id)) {
                        //f(x) = asin(x)
                        //f'(x) = 1/(2 sqrt(1-x^2)= 1/(pow((1-pow(x,2)),0.5)

                        ret->op_ = MULTIPLY;
                        ret->left_ = this->left_->Differentiate(id);


                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = DIVIDE;

                        ret->right_->left_ = new Expression<T > ();
                        ret->right_->left_->op_ = CONSTANT;
                        ret->right_->left_ ->value_ = T(1.0);

                        ret->right_->right_ = new Expression<T > ();
                        ret->right_->right_->op_ = POW;

                        ret->right_->right_->left_ = new Expression<T > ();
                        ret->right_->right_->left_->op_ = MINUS;

                        ret->right_->right_->left_->left_ = new Expression<T > ();
                        ret->right_->right_->left_->left_->op_ = CONSTANT;
                        ret->right_->right_->left_->left_->value_ = T(1.0);

                        ret->right_->right_->left_->right_ = new Expression<T > ();
                        ret->right_->right_->left_->right_->op_ = POW;
                        ret->right_->right_->left_->right_->left_ = this->left_->Clone();

                        ret->right_->right_->left_->right_->right_ = new Expression<T > ();
                        ret->right_->right_->left_->right_->right_->op_ = CONSTANT;
                        ret->right_->right_->left_->right_->right_->value_ = T(2.0);

                        ret->right_->right_->right_ = new Expression<T > ();
                        ret->right_->right_->right_->op_ = CONSTANT;
                        ret->right_->right_->right_->value_ = T(0.5);
                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case ACOS:

                    if (this->left_->HasID(id)) {
                        //f(x) = acos(x)
                        //f'(x) = -1/(sqrt(1-x^2) = -1/(pow((1-pow(x,2)),0.5)
                        //-1/sqrt(1-x^2)
                        ret->op_ = MULTIPLY;
                        ret->left_ = new Expression<T > ();
                        ret->left_->op_ = MULTIPLY;
                        ret->left_->left_ = new Expression<T > ();

                        ret->left_->left_->op_ = CONSTANT;
                        ret->left_->left_->value_ = T(-1.0);


                        ret->left_->right_ = this->left_->Differentiate(id);

                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = DIVIDE;

                        ret->right_->left_ = new Expression<T > ();
                        ret->right_->left_->op_ = CONSTANT;
                        ret->right_->left_ ->value_ = T(1.0);

                        ret->right_->right_ = new Expression<T > ();
                        ret->right_->right_->op_ = POW;

                        ret->right_->right_->left_ = new Expression<T > ();
                        ret->right_->right_->left_->op_ = MINUS;

                        ret->right_->right_->left_->left_ = new Expression<T > ();
                        ret->right_->right_->left_->left_->op_ = CONSTANT;
                        ret->right_->right_->left_->left_->value_ = T(1.0);

                        ret->right_->right_->left_->right_ = new Expression<T > ();
                        ret->right_->right_->left_->right_->op_ = POW;
                        ret->right_->right_->left_->right_->left_ = this->left_->Clone();

                        ret->right_->right_->left_->right_->right_ = new Expression<T > ();
                        ret->right_->right_->left_->right_->right_->op_ = CONSTANT;
                        ret->right_->right_->left_->right_->right_->value_ = T(2.0);

                        ret->right_->right_->right_ = new Expression<T > ();
                        ret->right_->right_->right_->op_ = CONSTANT;
                        ret->right_->right_->right_->value_ = T(0.5);
                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case ATAN:
                    if (this->left_->HasID(id)) {
                        //f(x) = atan(x)
                        //f'(x) 1/(x^2+1)

                        ret->op_ = DIVIDE;
                        ret->left_ = new Expression<T > ();
                        ret->left_->op_ = MULTIPLY;
                        ret->left_->right_ = new Expression<T > ();

                        ret->left_->right_->op_ = CONSTANT;
                        ret->left_->right_->value_ = T(1.0);


                        ret->left_->left_ = this->left_->Differentiate(id);

                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = PLUS;

                        ret->right_->left_ = new Expression<T > ();
                        ret->right_->left_->op_ = MULTIPLY;
                        ret->right_->left_->left_ = this->left_->Clone();
                        ret->right_->left_->right_ = this->left_->Clone();


                        ret->right_->right_ = new Expression<T > ();
                        ret->right_->right_->op_ = CONSTANT;
                        ret->right_->right_->value_ = T(1.0);

                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;

                    }
                case ATAN2:
                    //if w.r.t. check both expressions for id
                    if (this->left_->HasID(id)) {
                        //f(x) = atan2(x,y)
                        //f'(x) y/(x^2+y^2)

                        ret->op_ = DIVIDE;
                        ret->left_ = new Expression<T > ();
                        ret->left_->op_ = MULTIPLY;
                        ret->left_->left_ = this->right_->Clone(); //y
                        ret->left_->right_ = left_->Differentiate(id);


                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = PLUS;

                        ret->right_->left_ = new Expression<T > ();
                        ret->right_->left_->op_ = MULTIPLY;
                        ret->right_->left_->left_ = this->left_->Clone();
                        ret->right_->left_->right_ = this->left_->Clone();


                        ret->right_->right_ = new Expression<T > ();
                        ret->right_->right_->op_ = MULTIPLY;
                        ret->right_->right_->left_ = this->right_->Clone();
                        ret->right_->right_->right_ = this->right_->Clone();

                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case ATAN3:

                    //can be removed.
                    break;

                case ATAN4:
                    break;
                case SQRT:
                    if (this->left_->HasID(id)) {
                        //f(x) = sqrt(x)
                        //f'(x) = .5/sqrt(x)

                        ret->op_ = DIVIDE;
                        ret->left_ = new Expression<T > ();
                        ret->left_->op_ = MULTIPLY;

                        ret->left_->right_ = new Expression<T > ();
                        ret->left_->right_->value_ = T(0.5);

                        ret->left_->left_ = this->left_->Differentiate(id);


                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = SQRT;
                        ret->right_->left_ = this->left_->Clone();

                        //std::cout<<ret->ToString();
                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case POW:

                    if (this->left_->HasID(id)) {
                        //f(x) =  x^y
                        //f'(x) = yx^y-1

                        ret->op_ = MULTIPLY;

                        ret->left_ = new Expression<T > ();
                        ret->left_->op_ = MULTIPLY;
                        ret->left_->left_ = this->left_->Differentiate(id);
                        ret->left_->right_ = this->right_->Clone();



                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = POW;


                        ret->right_->left_ = this->left_->Clone();


                        ret->right_->right_ = new Expression<T > ();
                        ret->right_->right_->op_ = MINUS;
                        ret->right_->right_->left_ = this->right_->Clone();

                        ret->right_->right_->right_ = new Expression<T > ();
                        ret->right_->right_->right_->op_ = CONSTANT;
                        ret->right_->right_->right_->value_ = T(1.0);

                        //ret->Simplify();

                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
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
                    if (this->left_->HasID(id)) {
                        //f(x) = log(x)
                        //f'(x) = 1/x

                        ret->op_ = DIVIDE;
                        ret->left_ = new Expression<T > ();
                        ret->left_->op_ = MULTIPLY;
                        ret->left_->left_ = new Expression<T > ();
                        ret->left_->left_->op_ = CONSTANT;
                        ret->left_->left_->value_ = T(1.0);
                        ret->left_->right_ = this->left_->Differentiate(id);

                        ret->right_ = this->left_->Clone();


                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case LOG10:
                    //f(x) = log10(x)
                    //f'(x) = 1/(xlog(10))

                    if (this->left_->HasID(id)) {



                        ret->op_ = DIVIDE;

                        ret->left_ = new Expression<T > ();
                        ret->left_->op_ = MULTIPLY;

                        ret->left_->left_ = new Expression<T > ();
                        ret->left_->left_->op_ = CONSTANT;
                        ret->left_->left_->value_ = T(1.0);

                        ret->left_->right_ = this->left_->Differentiate(id);

                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = MULTIPLY;

                        ret->right_->left_ = this->left_->Clone();

                        ret->right_->right_ = new Expression<T > ();
                        ret->right_->right_->op_ = CONSTANT;
                        ret->right_->right_->value_ = log(T(10.0));
                        /*
                        ret->right_->right_->left_ = new Expression<T > ();
                        ret->right_->right_->left_->op_ = CONSTANT;
                        ret->right_->right_->left_->value_ = T(10.0);
                         */
                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = LOG;
                        ret->left_ = this->Clone();

                        //ret->Simplify();
                        return ret;
                    }
                case EXP:
                    //f(x) = e^x
                    //f'(x) =e^x

                    if (this->left_->HasID(id)) {

                        ret->op_ = MULTIPLY;
                        ret->left_ = this->left_->Differentiate(id);


                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = EXP;
                        ret->right_->left_ = this->left_->Clone();

                        //ret->Simplify();

                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case SINH:
                    if (this->left_->HasID(id)) {
                        //f(x) = sinh(x)
                        //f'(x) = cosh(x)

                        ret->op_ = MULTIPLY;
                        ret->left_ = this->left_->Differentiate(id);

                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = COSH;
                        ret->right_->left_ = this->left_->Clone();

                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case COSH:
                    if (this->left_->HasID(id)) {

                        ret->op_ = MULTIPLY;
                        ret->left_ = this->left_->Differentiate(id);

                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = SINH;
                        ret->right_->left_ = this->left_->Clone();

                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case TANH:
                    //f(x) = tanh(x)
                    //f'(x) =1- tanh(x)*tanh(x)


                    if (this->left_->HasID(id)) {

                        ret->op_ = MULTIPLY;

                        ret->left_ = this->left_->Differentiate(id);

                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = MULTIPLY;
                        ret->right_->left_ = new Expression<T > ();


                        ret->right_->left_->op_ = DIVIDE;
                        ret->right_->left_->left_ = new Expression<T > ();
                        ret->right_->left_->left_->op_ = CONSTANT;
                        ret->right_->left_->left_->value_ = T(1.0);


                        ret->right_->left_->right_ = new Expression<T > ();
                        ret->right_->left_->right_->op_ = COSH;
                        ret->right_->left_->right_->left_ = this->left_->Clone();


                        ret->right_->right_ = ret->right_->left_->Clone();
                        //
                        //                        ret->left_->op_ = MULTIPLY;
                        //                        ret->left_->right_ = new Expression<T > ();
                        //                        ret->left_->right_->op_ = CONSTANT;
                        //
                        //                        ret->left_->right_->value_ = T(1.0);
                        //
                        //                        ret->left_->left_ = this->left_->Differentiate(id);
                        //
                        //
                        //                        ret->right_ = new Expression<T > ();
                        //                        ret->right_->op_ = MULTIPLY;
                        //
                        //                        ret->right_->left_ = this->Clone();
                        //                        ret->right_->right_ = this->Clone();
                        //

                        //                ret->op_ = MULTIPLY;
                        //
                        //                ret->left_ = new Expression<T > ();
                        //                ret->left_->op_ = DIVIDE;
                        //                ret->left_->left_ = new Expression<T > ();
                        //                ret->left_->left_->op_ = CONSTANT;
                        //
                        //                ret->left_->right_ = new Expression<T > ();
                        //                ret->left_->right_->op_ = COSH;
                        //                ret->left_->right_->left_ = this->left_->Clone();
                        //
                        //
                        //                ret->right_ = new Expression<T > ();
                        //                ret->right_->op_ = DIVIDE;
                        //                ret->right_->left_ = new Expression<T > ();
                        //                ret->right_->left_->op_ = CONSTANT;
                        //
                        //                ret->right_->right_ = new Expression<T > ();
                        //                ret->right_->right_->op_ = COSH;
                        //                ret->right_->right_->left_ = this->left_->Clone();
                        //

                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }

                case FABS:

                    if (this->left_->HasID(id)) {

                        ret->op_ = DIVIDE;
                        ret->left_ = new Expression<T > ();
                        ret->left_->op_ = MULTIPLY;

                        ret->left_->left_ = this->left_->Differentiate(id);
                        ret->left_->right_ = this->left_->Clone();


                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = FABS;
                        ret->right_->left_ = this->left_->Clone();

                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
                        return ret;
                    }
                case FLOOR:
                    if (this->left_->id_ == id) {



                        ret->op_ = MULTIPLY;

                        ret->left_ = this->left_->Differentiate(id);

                        ret->right_ = new Expression<T > ();
                        ret->right_->op_ = FLOOR;
                        ret->right_->left_ = this->left_->Clone();

                        //ret->Simplify();
                        return ret;
                    } else {
                        ret->op_ = CONSTANT;
                        ret->value_ = T(0.0);

                        //ret->Simplify();
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
        void Solve(const uint32 &id){
            
        }
        
        /*!
         * Simplifies this expression, effectively reducing
         * the amount of memory needed to store the expression.
         */
        void Simplify() {
            //            // std::cout<<"SIMPLIFY\n\n";
            //            Expression<T> * temp;
            //
            //
            //            if (this->left_) {
            //                this->left_->Simplify();
            //            }
            //
            //            if (this->right_) {
            //                this->right_->Simplify();
            //            }
            //
            //
            //            if (this->left_ != NULL && this->right_ != NULL) {
            //
            //
            //
            //
            //
            //                switch (this->op_) {
            //
            //                    case MINUS:
            //                        if (this->left_->op_ == CONSTANT && this->right_->op_ == CONSTANT) {
            //
            //                            this->op_ = CONSTANT;
            //
            //                            this->value_ = this->left_->value_ - this->right_->value_;
            //
            //                            delete this->left_;
            //                            delete this->right_;
            //                        }
            //
            //
            //
            //                        break;
            //
            //                    case PLUS:
            //                        if (this->left_->op_ == CONSTANT && this->left_->value_ == T(0)) {
            //
            //                            temp = this->right_;
            //                            delete this->left_;
            //                            this->op_ = this->right_->op_;
            //                            this->id_ = this->right_->id_;
            //                            this->value_ = this->right_->value_;
            //
            //                            this->left_ = temp->left_->Clone();
            //                            this->right_ = temp->right_->Clone();
            //                            delete temp;
            //
            //                        } else if (this->right_->op_ == CONSTANT && this->right_->value_ == T(0)) {
            //
            //                            temp = this->left_;
            //                            delete this->right_;
            //                            this->op_ = this->left_->op_;
            //                            this->id_ = this->left_->id_;
            //                            this->value_ = this->left_->value_;
            //
            //                            this->left_ = temp->left_->Clone();
            //                            this->right_ = temp->right_->Clone();
            //                            delete temp;
            //
            //
            //                        } else if (this->left_->op_ == CONSTANT && this->right_->op_ == CONSTANT) {
            //
            //                            this->op_ = CONSTANT;
            //
            //                            this->value_ = this->left_->value_ + this->right_->value_;
            //
            //                            delete this->left_;
            //                            delete this->right_;
            //                        }
            //
            //
            //                        break;
            //
            //                    case DIVIDE:
            //
            //                        if (this->left_->op_ == CONSTANT && this->right_->op_ == CONSTANT) {
            //
            //                            this->op_ = CONSTANT;
            //
            //                            this->value_ = this->left_->value_ / this->right_->value_;
            //
            //                            delete this->left_;
            //                            delete this->right_;
            //                        }
            //
            //
            //                        break;
            //
            //                    case MULTIPLY:
            //
            //                        //                                if (this->left_->op_ == CONSTANT && this->left_->value_ == T(0.0)) {
            //                        //
            //                        //
            //                        //                                    this->op_ = CONSTANT;
            //                        //                                
            //                        //                                    this->value_ = T(0);
            //                        //
            //                        //
            //                        //                                    delete this->left_;
            //                        //                                    delete this->right_;
            //                        //
            //                        //                                } else if (this->right_->op_ == CONSTANT && this->right_->value_ == T(0.0)) {
            //                        //
            //                        //
            //                        //                                    this->op_ = CONSTANT;
            //                        //                                    
            //                        //                                    this->value_ = T(0);
            //                        //
            //                        //
            //                        //                                    delete this->left_;
            //                        //                                    delete this->right_;
            //                        //
            //                        //                                } else if (this->left_->op_ == CONSTANT && this->left_->value_ == T(1.0)) {
            //                        //
            //                        //
            //                        //                                    temp = this->right_;
            //                        //
            //                        //                                    this->op_ = temp->op_;
            //                        //                                    this->value_ = temp->value_;
            //                        //                                   
            //                        //                                    this->id_ = temp->id_;
            //                        //                                    this->left_ = temp->left_->Clone();
            //                        //                                    this->right_ = temp->right_->Clone();
            //                        //
            //                        //
            //                        //                                    delete this->left_;
            //                        //                                    delete temp;
            //                        //
            //                        //
            //                        //                                } else if (this->right_->op_ == CONSTANT && this->right_->value_ == T(1.0)) {
            //                        //
            //                        //
            //                        //                                    temp = this->left_;
            //                        //
            //                        //                                    this->op_ = temp->op_;
            //                        //                                    this->value_ = temp->value_;
            //                        //                                    
            //                        //                                    this->id_ = temp->id_;
            //                        //                                    this->left_ = temp->left_->Clone();
            //                        //                                    this->right_ = temp->right_->Clone();
            //                        //
            //                        //
            //                        //                                    delete this->right_;
            //                        //                                    delete temp;
            //                        //
            //                        //
            //                        //                                } else if (this->left_->op_ == CONSTANT && this->right_->op_ == CONSTANT) {
            //                        //
            //                        //                                    this->op_ = CONSTANT;
            //                        //                                 
            //                        //
            //                        //                                    this->value_ = this->left_->value_ * this->right_->value_;
            //                        //
            //                        //                                    delete this->left_;
            //                        //                                    delete this->right_;
            //                        //                                }
            //                        //
            //                        //                                break;
            //                    case SIN:
            //
            //
            //                        break;
            //                    case COS:
            //
            //                        break;
            //                    case TAN:
            //                        break;
            //                    case ASIN:
            //                        break;
            //                    case ACOS:
            //                        break;
            //                    case ATAN:
            //                        break;
            //                    case ATAN2:
            //                        break;
            //                    case ATAN3:
            //                        break;
            //                    case ATAN4:
            //                    case SQRT:
            //                        break;
            //                    case POW:
            //
            //                        if (right_->Evaluate() == T(1)) {
            //                            std::cout << "eguals 1....\n";
            //                        }
            //
            //                        if (right_->Evaluate() == T(0)) {
            //                            std::cout << "eguals 0....\n";
            //                        }
            //
            //
            //                        break;
            //                    case POW1:
            //                        break;
            //                    case POW2:
            //                        break;
            //                    case LOG:
            //                        break;
            //                    case LOG10:
            //                        break;
            //                    case EXP:
            //                        break;
            //                    case SINH:
            //                        break;
            //                    case COSH:
            //                        break;
            //                    case TANH:
            //                        break;
            //                    case NONE:
            //                        break;
            //                    default:
            //                        break;
            //                }
            //            }

        }

        /*!
         * Return a list of differentiable ids in this
         * expression.
         * 
         * @param vars
         */
        void VariableIds(std::vector< uint32_t> &vars) {

            if (this->left_ != NULL) {
                this->left_->VariableIds(vars);
            }


            if (this->right_ != NULL) {
                this->right_->VariableIds(vars);
            }

            if (this->op_ == VARIABLE) {
                bool exists = false;
                for (size_t i = 0; i < vars.size(); i++) {
                    if (vars.at(i) == this->id_) {
                        exists = true;
                    }
                }

                if (!exists) {
                    vars.push_back(this->id_);
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

            if (this->left_ != NULL) {
                l = this->left_->ToString(vars);
            }
            temp.str("");

            if (this->right_ != NULL) {
                r = this->right_->ToString(vars);
            }
            ss << "(";

            switch (this->op_) {
                case CONSTANT:
                    ss << this->value_ << "";
                    break;
                case VARIABLE:
                    ss << "x" << this->id_;


                    temps << "x" << this->id_ << "/*" << this->value_ << "*/";

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
        std::string ToPrettyString() {
            std::stringstream ss;
            std::stringstream temp;

            std::string l, r;


            if (this->left_ != NULL) {
                l = this->left_->ToPrettyString();
            }
            temp.str("");

            if (this->right_ != NULL) {
                r = this->right_->ToPrettyString();
            }
            ss << "";

            switch (this->op_) {
                case CONSTANT:
                    ss << "CONST[" << this->value_ << "]";
                    break;
                case VARIABLE:
                    ss << "VAR[" << this->value_ << ",ID[" << this->id_ << "]" << "]";
                    break;
                case MINUS:
                    ss << "" << l << "\n - \n" << r;
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
                case NONE:
                    break;
                default:
                    break;
            }

            ss << "";

            return ss.str();

        }

        /*!
         * Represent this expression as a string. ADNumbers are represented
         * in wkt format by value and id. Constants are represented by value.
         *
         */
        std::string ToString() {
            std::stringstream ss;
            std::stringstream temp;

            std::string l, r;


            if (this->left_ != NULL) {
                l = this->left_->ToString();
            }
            temp.str("");

            if (this->right_ != NULL) {
                r = this->right_->ToString();
            }
            ss << "(";

            switch (this->op_) {
                case CONSTANT:
                    ss << "CONST[" << this->value_ << "]";
                    break;
                case VARIABLE:
                    ss << "VAR[" << this->value_ << ",ID[" << this->id_ << "]" << "]";
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

            return id_;
        }

        /*!
         * Set the unique identifier associated with this expression.
         * @param id
         */
        void SetId(unsigned long id) {

            this->id_ = id;
        }

        /*!
         * Returns the left branch of this expression tree.
         * @return 
         */
        Expression<T> GetLeft() const {

            return left_;
        }

        /*!
         * Sets the left branch of this expression tree.
         *
         */
        void SetLeft(Expression<T> *left) {

            this->left_ = left;
        }

        /*!
         * Return the operator for this expression.
         * @return 
         */
        Operation GetOp() const {

            return op_;
        }

        /*!
         * Set the operator for this expression.
         * 
         */
        void SetOp(const Operation &op) {

            this->op_ = op;
        }

        /*!
         * Returns the right branch of this expression tree.
         * @return 
         */
        Expression<T> GetRight() const {

            return right_;
        }

        /*!
         * Sets the right branch of this expression tree.
         * 
         */
        void SetRight(Expression<T> *right) {

            this->right_ = right;
        }

        /*!
         * Returns the value assigned to this 
         * expression.
         * @return 
         */
        T GetValue() const {

            return value_;
        }

        /*!
         * Sets the value assigned to this 
         * expression.
         * 
         */
        void SetValue(T value) {
            this->value_ = value;
        }



    private:
        T epsilon_m;
        Expression<T>* left_; //left branch
        Expression<T>* right_; //right branch
        T value_; //raw value
        unsigned int id_; //unique id
        Operation op_; //operation





    };

    /*!
     * Template class ADNumber. Automatic differentiation class.
     * Supports forward and reverse mode traversal of the chain rule.
     * Supports higher order derivatives, as well as the partial nth
     * derivative with respect to a ADNumber.
     * 
     */
    template<class T>
    class ADNumber {
    public:

        /*!
         * Default constructor.
         */
        ADNumber() :
        expression_(new Expression<T>()),
        value_(T(0.0)),
        fderivative_(T(1)),
        variableName_(std::string("na")),
        id_(IDGenerator::instance()->next()) {
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
        ADNumber(T value, T derivative = T(1.0)) :
        expression_(new Expression<T>()),
        value_(value),
        fderivative_(derivative),
        variableName_(std::string("x")),
        id_(IDGenerator::instance()->next()) {
            this->Initialize();

        }

        /*!
         * Constructor
         * 
         * @param name
         * @param value
         * @param derivative
         */
        ADNumber(std::string name, T value, T derivative = T(1.0)) :
        value_(value),
        fderivative_(derivative),
        variableName_(name),
        expression_(new Expression<T>()),
        id_(IDGenerator::instance()->next()) {
            this->Initialize();
        }

        /*!
         * Copy Constructor.
         * 
         * @param orig
         */
        ADNumber(const ADNumber& orig) :
        value_(T(0.0)),
        fderivative_(T(1)),
        variableName_(std::string("na")),
        id_(orig.GetID()) {

            this->variableName_ = orig.variableName_;
            this->value_ = orig.value_;
            this->fderivative_ = orig.fderivative_;
            this->id_ = orig.id_;

            if (orig.expression_ != NULL) {

                this->expression_ = orig.expression_->Clone();
            }
            // this->Initialize();
        }

        //        /*!
        //         * Copy Constructor.
        //         * 
        //         * @param orig
        //         */
        //        ADNumber(const ADNumber* orig) :
        //        value_(T(0)),
        //        fderivative_(T(1)),
        //        variableName_(std::string("na")),
        //        id_(orig->GetID()) {
        //
        //            std::cout << "copy * called!!!!\n";
        //
        //            this->variableName_ = orig->variableName_;
        //            this->value_ = orig->value_;
        //            this->fderivative_ = orig->fderivative_;
        //            this->id_ = orig->GetID();
        //
        //            if (orig->expression_ != NULL) {
        //                this->expression_ = orig->expression_->Clone();
        //            }
        //            //this->Initialize();
        //        }

#ifdef USE_MEMORY_POOL

        void* operator new (size_t size) throw (std::bad_alloc) {
            assert(size == sizeof (ADNumber));
            void* ptr = MemoryPool < sizeof (ADNumber)>::instance().allocate();
            return ptr;
        }

        void operator delete (void* ptr)throw () {
            // void *ptr = MemoryPool<sizeof (Test)>::instance().allocate();
            MemoryPool < sizeof (ADNumber)>::instance().deallocate((void*) ptr);
        }
#endif

        /*!
         * Destructor.
         */
        virtual ~ADNumber() {

            if (this->expression_ != NULL) {
                delete this->expression_;
            }
        }

        /*!
         * returns this value.
         */
        operator T() const {
            return this->value_;
        }

        /*!
         * Returns this value.
         */
        operator ADNumber<T>() const {
            return this;
        }

        /*!
         * In member assignment operator to set this 
         * equal to another ADNumber val.
         */
        ADNumber<T> & operator =(const ADNumber<T> &val) {
            if (val.expression_ != NULL) {
                delete this->expression_;
            }

            this->expression_ = val.expression_->Clone();
            this->value_ = val.GetValue();
            this->fderivative_ = val.Forward();
            this->id_ = val.id_;
            this->variableName_ = val.variableName_;

            return *this;
        }

        /*!
         * In member assignment operator to set this value 
         * equal to val with derivative set to 1.
         */
        ADNumber<T> & operator =(const T & val) {
            this->value_ = val;
            this->fderivative_ = T(1.0);
            this->id_ = uint32_t(IDGenerator::instance()->next());
            delete this->expression_;
            this->expression_ = new Expression<T > ();


            this->Initialize();

            return *this;
        }

        /*!
         * In member addition operator.
         * Returns ADNumber<T> with this value + rhs value &
         * this derivative + rhs derivitive.
         */
        ADNumber<T> operator +(const ADNumber<T>& rhs) const {
            ADNumber<T > ret(T(this->value_ + rhs.GetValue()),
                    this->fderivative_ + rhs.Forward());

            ret.expression_->SetOp(PLUS);
            ret.expression_->SetLeft(this->expression_->Clone());
            ret.expression_->SetRight(rhs.expression_->Clone());

            return ret;
        }

        /*!
         * In member addition operator.
         * Returns ADNumber<T> with this value + rhs value &
         * this derivative + 0.
         */
        ADNumber<T> operator +(const T & rhs) const {
            ADNumber<T > ret(T(this->value_ + rhs),
                    T(this->fderivative_));

            ret.expression_->SetOp(PLUS);
            ret.expression_->SetLeft(this->expression_->Clone());

            Expression<T> *temp = new Expression<T > ();
            temp->SetOp(CONSTANT);
            temp->SetValue(rhs);
            ret.expression_->SetRight(temp);

            return ret;
        }

        /*!
         * In member subtraction operator.
         * Returns ADNumber<T> with this value - rhs value &
         * this derivative - rhs derivitive.
         *
         */
        ADNumber<T> operator -(const ADNumber<T>& rhs) const {
            ADNumber<T > ret(T(this->value_ - rhs.GetValue()),
                    T(this->fderivative_ - rhs.Forward()));


            ret.expression_->SetOp(MINUS);
            ret.expression_->SetLeft(this->expression_->Clone());
            ret.expression_->SetRight(rhs.expression_->Clone());

            return ret;
        }

        /*!
         * In member subtraction operator.
         * Returns ADNumber<T>(this value - rhs value, this derivative - 0).
         */
        ADNumber<T> operator -(const T & rhs) const {
            ADNumber<T > ret(T(this->value_ - rhs),
                    T(this->fderivative_));

            ret.expression_->SetOp(MINUS);
            ret.expression_->SetLeft(this->expression_->Clone());

            Expression<T> *temp = new Expression<T > ();
            temp->SetOp(CONSTANT);
            temp->SetValue(rhs);
            ret.expression_->SetRight(temp);

            return ret;
        }

        /*!
         * In member multiplication operator.
         * Returns ADNumber<T>(this value * rhs value,
         * this value_ * rhs derivative  + rhs value * this derivative).
         */
        ADNumber<T> operator *(const ADNumber<T>& rhs) const {
            ADNumber<T > ret(T(this->value_ * rhs.GetValue()),
                    T(this->value_ * rhs.Forward() +
                    rhs.GetValue() * this->fderivative_));

            ret.expression_->SetOp(MULTIPLY);
            ret.expression_->SetLeft(this->expression_->Clone());
            ret.expression_->SetRight(rhs.expression_->Clone());

            return ret;
        }

        /*!
         * In member multiplication operator.
         * Returns ADNumber<T>(this value * rhs value,
         * this value_ * 0  + rhs value * this derivative).
         */
        ADNumber<T> operator *(const T & rhs) const {
            ADNumber<T > ret(T(this->value_ * rhs),
                    T(this->value_ * 0 + rhs * this->fderivative_));

            ret.expression_->SetOp(MULTIPLY);
            ret.expression_->SetLeft(this->expression_->Clone());

            Expression<T> *temp = new Expression<T > ();
            temp->SetOp(CONSTANT);
            temp->SetValue(rhs);
            ret.expression_->SetRight(temp);

            return ret;
        }

        /*!
         * In member division operator.
         * Returns ADNumber<T>(this value * rhs value,
         * (this value_ * rhs derivative  - rhs value * this derivative)/(rhs value * rhs value)).
         */
        ADNumber<T> operator /(const ADNumber<T>& rhs) const {
            ADNumber<T > ret(T(this->value_ / rhs.value_),
                    T((rhs.GetValue() * this->fderivative_ -
                    this->value_ * rhs.Forward()) / (rhs.GetValue() * rhs.GetValue())));

            ret.expression_->SetOp(DIVIDE);
            ret.expression_->SetLeft(this->expression_->Clone());
            ret.expression_->SetRight(rhs.expression_->Clone());

            return ret;
        }

        /*!
         * In member division operator.
         * Returns ADNumber<T>(this value * rhs value,
         * (this value_ * rhs derivative  - rhs value * 0)/(rhs value * rhs value)).
         */
        ADNumber<T> operator /(const T & rhs) const {
            ADNumber<T > ret(T(this->GetValue() / rhs),
                    T((rhs * this->GetValue() - this->GetValue() * 0) / (rhs * rhs)));

            ret.expression_->SetOp(DIVIDE);
            ret.expression_->SetLeft(this->expression_->Clone());

            Expression<T> *temp = new Expression<T > ();
            temp->SetOp(CONSTANT);
            temp->SetValue(rhs);
            ret.expression_->SetRight(temp);

            return ret;
        }

        /*!
         * In member addition assignment operator.
         * @param rhs
         * @return 
         */
        ADNumber<T> operator +=(const ADNumber<T>& rhs) {
            ADNumber<T> ret = (*this +rhs);
            if (ret.expression_ != NULL) {
                delete this->expression_;
                this->expression_ = ret.expression_->Clone();
            }
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();


            return ret;
        }

        /*!
         * In member addition subtraction operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T> operator -=(const ADNumber<T>& rhs) {
            ADNumber<T> ret = (*this -rhs);
            if (ret.expression_ != NULL) {
                delete this->expression_;
                this->expression_ = ret.expression_->Clone();
            }
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();

            return ret;
        }

        /*!
         * In member multiplication assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T> operator *=(const ADNumber<T>& rhs) {
            ADNumber<T> ret = (*this * rhs);
            if (ret.expression_ != NULL) {
                delete this->expression_;
                this->expression_ = ret.expression_->Clone();
            }
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();

            return ret;
        }

        /*!
         * In member division assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T> operator /=(const ADNumber<T>&rhs) {
            ADNumber<T> ret = (*this / rhs);
            if (ret.expression_ != NULL) {
                delete this->expression_;
                this->expression_ = ret.expression_->Clone();
            }
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();

            return ret;
        }

        /*!
         * In member addition assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T> operator +=(const T & rhs) {
            ADNumber<T> ret = (*this +rhs);
            if (ret.expression_ != NULL) {
                delete this->expression_;
                this->expression_ = ret.expression_->Clone();
            }
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();

            return ret;
        }

        /*!
         * In member subtraction assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T> operator -=(const T & rhs) {
            ADNumber<T> ret = (*this -rhs);
            if (ret.expression_ != NULL) {
                delete this->expression_;
                this->expression_ = ret.expression_->Clone();
            }
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();

            return ret;
        }

        /*!
         * In member multiplication assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T> operator *=(const T & rhs) {
            ADNumber<T> ret = (*this * rhs);
            if (ret.expression_ != NULL) {
                delete this->expression_;
                this->expression_ = ret.expression_->Clone();
            }
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();

            return ret;
        }

        /*!
         * In member division assignment operator.
         * 
         * @param rhs
         * @return 
         */
        ADNumber<T> operator /=(const T & rhs) {
            ADNumber<T> ret = (*this / rhs);
            if (ret.expression_ != NULL) {
                delete this->expression_;
                this->expression_ = ret.expression_->Clone();
            }
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();

            return ret;
        }

        /*!
         * In member suffix increment operator.
         * 
         * @return 
         */
        ADNumber<T> operator ++() {

            ADNumber<T> ret(*this+T(1.0));
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();
            delete this->expression_;
            this->expression_ = ret.expression_->Clone();

            return *this;
        }

        /*!
         * In member suffix decrement operator.
         * 
         * @return 
         */
        ADNumber<T> operator --() {
            ADNumber<T> ret(*this-T(1.0));
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();
            delete this->expression_;
            this->expression_ = ret.expression_->Clone();

            return *this;
        }

        /*!
         * In member prefix increment operator.
         * 
         * @param 
         */
        ADNumber<T> operator ++(int) {
            ADNumber<T> temp = *this;

            ADNumber<T> ret(*this+T(1.0));
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();
            delete this->expression_;
            this->expression_ = ret.expression_->Clone();

            return temp;

        }

        /*!
         * In member prefix decrement operator.
         * 
         * @param 
         */
        ADNumber<T> operator --(int) {
            ADNumber<T> temp = *this;

            ADNumber<T> ret(*this-T(1.0));
            this->value_ = ret.GetValue();
            this->fderivative_ = ret.Forward();
            delete this->expression_;
            this->expression_ = ret.expression_->Clone();

            return *temp;
        }

        /*!
         * Returns the computed value.
         */
        const T GetValue() const {

            return this->value_;
        }

        /*!
         * Returns the result of all derivatives in the expression chain via forward
         * accumulation. Forward mode automatic differentiation is accomplished by
         * the dual numbers method.
         * Source:http://en.wikipedia.org/wiki/Automatic_differentiation#Automatic_differentiation_using_dual_numbers
         */
        const T Forward() const {
            return this->fderivative_;
        }

        /*!
         * Returns the result of all derivatives in the expression chain via reverse
         * accumulation.  Reverse accumulation traverses the chain rule from left
         * to right, or in the case of the computational graph,
         * from top to bottom.
         * Source: http://en.wikipedia.org/wiki/Automatic_differentiation#Reverse_accumulation
         */
        const T Reverse() const {
            if (this->expression_ == NULL) {
                return T(0);
            }

            Expression<T> * exp = this->expression_->Differentiate();
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
            // Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            if (this->expression_ == NULL) {
                return T(0);
            } else {
                if (vars.size() == 0) {
                    return this->GetValue();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_->Differentiate(vars.at(0)->GetID());

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
            // Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            Expression<T> *temp;
            Expression<T> *exp =
                    this->expression_->Differentiate(var0.GetID());

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
            if (this->expression_ == NULL) {
                return T(0);
            } else {
                if (vars.size() == 0) {
                    return this->GetValue();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_->Differentiate(vars.at(0)->GetID());

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


            if (this->expression_ == NULL) {
                return T(0);
            } else {
                if (order == 0) {
                    return this->GetValue();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_->Differentiate();

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


            if (this->expression_ == NULL) {
                return T(0);
            } else {
                if (order == 0) {
                    return this->expression_->PropagatedError();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_->Differentiate();

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

            if (!this->expression_) {
                return T(0);
            } else {
                if (order == 0) {
                    return this->GetValue();
                }


                Expression<T> *temp;
                Expression<T> *exp = this->expression_->Differentiate(wrt.GetID());


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

            if (this->expression_ == NULL) {
                return T(0);
            } else {
                if (order == 0) {
                    return this->NthError(0);
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_->Differentiate(wrt.GetID());

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

            std::vector<uint32_t> vars;
            this->expression_->VariableIds(vars);


            T temp = T(0);
            T squared_epsilon = std::numeric_limits<T>::epsilon() * std::numeric_limits<T>::epsilon();
            T dif;
            for (size_t i = 0; i < vars.size(); i++) {
                Expression<T> *exp =
                        this->expression_->Differentiate(vars.at(i));
                dif = exp->Evaluate();
                temp += dif * dif * squared_epsilon;
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

            if (this->expression_ == NULL) {
                return "NA";
            } else {
                if (order == 0) {
                    return this->expression_->ToString();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_->Differentiate();

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

            // std::cout << __func__ << ": Not yet Implemented!\n";
            //return "Not yet Implemented!\n";
            std::vector<std::string> args;

            std::string ret; // = this->expression_->ToString(args);


            if (this->expression_ == NULL) {
                return "NA";
            } else {
                if (order == 0) {
                    ret = this->expression_->ToString(args);
                } else {

                    Expression<T> *temp;
                    Expression<T> *exp = this->expression_->Differentiate();

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
         * Return the nth order partial derivative as a std::string.
         */
        const std::string NthPartialToCPPFunction(std::string name, const ADNumber<T> &wrt, const unsigned int &order) {

            // std::cout << __func__ << ": Not yet Implemented!\n";
            //return "Not yet Implemented!\n";
            std::vector<std::string> args;

            std::string ret; // = this->expression_->ToString(args);

            if (this->expression_ == NULL) {
                return "NA";
            } else {
                if (order == 0) {
                    ret = this->expression_->ToString(args);
                } else {

                    Expression<T> *temp;
                    Expression<T> *exp = this->expression_->Differentiate(wrt.GetID());

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

            // std::cout << __func__ << ": Not yet Implemented!\n";
            //return "Not yet Implemented!\n";

            if (this->expression_ == NULL) {
                return "NA";
            } else {
                if (order == 0) {
                    return this->expression_->ToString();
                }

                Expression<T> *temp;
                Expression<T> *exp = this->expression_->Differentiate(wrt.GetID());

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

            return this->variableName_;
        }

        void SetName(const std::string &name) {
            this->variableName_ = name;
        }

        /*
         * Return the unique identifier for this ADNumber.
         */
        const uint32_t GetID() const {

            return this->id_;
        }


        //expression tree, used for reverse mode calculation.
        Expression<T> *expression_;

    protected:
        //computed value.
        T value_;

        //forward computed derivative.(direct)
        T fderivative_;


    private:
        //initialize this expression.

        void Initialize() {

            this->expression_->SetValue(this->value_);
            this->expression_->SetId(this->id_);
            this->expression_->SetOp(VARIABLE);


        }


        //reverse computed derivative.holds the expression tree calculation to
        //avoid repeated evaluations of the expression tree.(adjoint)
        //T rderivative_;

        //reverse flag
        bool reverse_computed;

        //Default is x
        std::string variableName_;

        //unique id
        uint32_t id_;



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

        ret.expression_->SetOp(MINUS);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(PLUS);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(DIVIDE);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(MULTIPLY);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(MINUS);
        ret.expression_->SetLeft(exp);
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(PLUS);
        ret.expression_->SetLeft(exp);
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(DIVIDE);
        ret.expression_->SetLeft(exp);
        ret.expression_->SetRight(rhs.expression_->Clone());

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
                T(lhs * rhs.Forward() + rhs.GetValue() * 0));

        Expression<T> *exp = new Expression<T > ();
        exp->SetValue(lhs);
        exp->SetOp(CONSTANT);

        ret.expression_->SetOp(MULTIPLY);
        ret.expression_->SetLeft(exp);
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(MINUS);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(exp);

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

        ret.expression_->SetOp(PLUS);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(exp);

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

        ret.expression_->SetOp(DIVIDE);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(exp);

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

        ret.expression_->SetOp(MULTIPLY);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(exp);

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

        ret.expression_->SetOp(ad::ATAN);
        ret.expression_->SetLeft(val.expression_->Clone());

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

        ret.expression_->SetOp(ad::ATAN2);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(ad::ATAN2);
        ret.expression_->SetLeft(exp);
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(ad::ATAN2);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(exp);

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

        ret.expression_->SetOp(ad::COS);
        ret.expression_->SetLeft(val.expression_->Clone());

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

        ret.expression_->SetOp(ad::EXP);
        ret.expression_->SetLeft(val.expression_->Clone());

        return ret;
    }

    /*!
     * Compute natural logarithm of val.
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> log(const ad::ADNumber<T> &val) {
        ad::ADNumber<T> ret(log(val.GetValue()), T(1.0) / val.GetValue());

        ret.expression_->SetOp(ad::LOG);
        ret.expression_->SetLeft(val.expression_->Clone());

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



        ret.expression_->SetOp(ad::LOG10);
        ret.expression_->SetLeft(val.expression_->Clone());


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

        ret.expression_->SetOp(ad::POW);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(ad::POW);
        ret.expression_->SetLeft(exp);
        ret.expression_->SetRight(rhs.expression_->Clone());

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

        ret.expression_->SetOp(ad::POW);
        ret.expression_->SetLeft(lhs.expression_->Clone());
        ret.expression_->SetRight(exp);

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

        ret.expression_->SetOp(ad::SIN);
        ret.expression_->SetLeft(val.expression_->Clone());

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
        ad::Expression<T>* right = new ad::Expression<T > ();
        right->SetOp(ad::CONSTANT);
        right->SetValue(T(0.5));

        //just use pow!!!
        ret.expression_->SetOp(ad::POW);
        ret.expression_->SetLeft(val.expression_->Clone());
        ret.expression_->SetRight(right);
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

        ret.expression_->SetOp(ad::TAN);
        ret.expression_->SetLeft(val.expression_->Clone());

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

        ret.expression_->SetOp(ad::ACOS);
        ret.expression_->SetLeft(val.expression_->Clone());

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

        ret.expression_->SetOp(ad::ASIN);
        ret.expression_->SetLeft(val.expression_->Clone());

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

        ret.expression_->SetOp(ad::SINH);
        ret.expression_->SetLeft(val.expression_->Clone());
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

        ret.expression_->SetOp(ad::COSH);
        ret.expression_->SetLeft(val.expression_->Clone());

        return ret;
    }

    /*!
     * Returns the hyperbolic tangent of val.
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> tanh(const ad::ADNumber<T> &val) {
        T temp = cosh(val.GetValue());
        ad::ADNumber<T> ret(tanh(val.GetValue()));

        ret.expression_->SetOp(ad::TANH);
        ret.expression_->SetLeft(val.expression_->Clone());

        return ret;
    }

    /*!
     * Compute absolute value of val.
     * @param val
     * @return 
     */
    template<class T> ad::ADNumber<T> fabs(const ad::ADNumber<T> &val) {

        ad::ADNumber<T> ret(fabs(val.GetValue()), fabs(val.Forward()));

        ret.expression_->SetOp(ad::FABS);
        ret.expression_->SetLeft(val.expression_->Clone());

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
        ret.expression_->SetOp(ad::FLOOR);
        ret.expression_->SetLeft(val.expression_->Clone());
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
        std::cout<<"solve not yet implemented....\n";
         ad::ADNumber<T> ret;
         
         if(lhs.expression_->HasID(var.GetID())){
             std::cout<<"left side contains var...\n";
         }
         
         if(rhs.expression_->HasID(var.GetID())){
             std::cout<<"left side contains var...\n";
         }
         
         
         return ret;

    }


    //namespace std {

    template<class T>
    class numeric_limits<ad::ADNumber<T> > : public numeric_limits<T> {
    };





}


template<class T> 
inline T diff(const ad::ADNumber<T> &var, unsigned int order = 1){
    return var.Nth(order);
}


template<class T> 
inline T diff(const ad::ADNumber<T> &var, const ad::ADNumber<T> &wrt, unsigned int order = 1){
    return var.NthPartial(wrt,order);
}

typedef ad::ADNumber<double> dvar;
typedef ad::ADNumber<float> fvar;

#endif	/* ADNUMBER_HPP */
