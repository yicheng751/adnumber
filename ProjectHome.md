Software to compute derivatives of arbitrary order to machine precision. Support for nth partial and mixed derivatives.



# Examples #

## Derivatives ##

```

#include <iostream>
#include "ADNumber.hpp"

using namespace std;

typedef ad::ADNumber<double> var;

/*
 * Example using ADNumber.
 * 
 * reference: http://www.wolframalpha.com/input/?i=derivative+of+x%5E4+sin+x&lk=3
 */
int main(int argc, char** argv) {

    
    var f;
    var x = 1;
    f = std::pow(x, 4.0) * sin(x);
    
    std::cout << "d/dx(" << f.NthToString(0)
            << ") = " << f.NthToString(1)
            << " = " << f.Forward() << std::endl;


    return 0;
}

```

**Output**


```
d/dx(((pow((VAR[1,ID[2]]),(CONST[4]))) * (sin((VAR[1,ID[2]]))))) = (((((CONST[1]) * (CONST[4])) * (pow((VAR[1,ID[2]]),((CONST[4]) - (CONST[1]))))) * (sin((VAR[1,ID[2]])))) + ((pow((VAR[1,ID[2]]),(CONST[4]))) * (cos((VAR[1,ID[2]]))))) = 3.90619
```

**or**


```

#include <iostream>
#include "ADNumber.hpp"

using namespace std;

typedef ad::ADNumber<double> var;

/*
 * Example using ADNumber.
 * 
 * reference: http://www.wolframalpha.com/input/?i=derivative+of+x%5E4+sin+x&lk=3
 */
int main(int argc, char** argv) {

    
    var f;
    var x = 1;
    f = std::pow(x, 4.0) * sin(x);
    
    std::cout << "d/dx(" << f.NthToString(0)
            << ") = " << f.NthToString(1)
            << " = " << f.Reverse()<< std::endl;


    return 0;
}

```

**Output**


```
d/dx(((pow((VAR[1,ID[2]]),(CONST[4]))) * (sin((VAR[1,ID[2]]))))) = (((((CONST[1]) * (CONST[4])) * (pow((VAR[1,ID[2]]),((CONST[4]) - (CONST[1]))))) * (sin((VAR[1,ID[2]])))) + ((pow((VAR[1,ID[2]]),(CONST[4]))) * (cos((VAR[1,ID[2]]))))) = 3.90619
```

**or**

```
#include <iostream>
#include "ADNumber.hpp"

using namespace std;

typedef ad::ADNumber<double> var;

/*
 * Example using ADNumber.
 * 
 * reference: http://www.wolframalpha.com/input/?i=derivative+of+x%5E4+sin+x&lk=3
 */
int main(int argc, char** argv) {

    
    var f;
    var x = 1;
    f = std::pow(x, 4.0) * sin(x);
    
    std::cout << "d/dx(" << f.NthToString(0)
            << ") = " << f.NthToString(1)
            << " = " << f.Nth(0)<< std::endl;


    return 0;
}
```

**Output**


```
d/dx(((pow((VAR[1,ID[2]]),(CONST[4]))) * (sin((VAR[1,ID[2]]))))) = (((((CONST[1]) * (CONST[4])) * (pow((VAR[1,ID[2]]),((CONST[4]) - (CONST[1]))))) * (sin((VAR[1,ID[2]])))) + ((pow((VAR[1,ID[2]]),(CONST[4]))) * (cos((VAR[1,ID[2]]))))) = 3.90619
```



## Higher Derivatives ##
```

#include <iostream>
#include "ADNumber.hpp"

using namespace std;

typedef ad::ADNumber<double> var;

/*
 * Example of higher derivatives.
 * 
 * reference: http://www.wolframalpha.com/input/?i=second+derivative+of+sin%282x%29&lk=3
 */
int main(int argc, char** argv) {

    
    var f;
    var x = 1;
   
    
    f =  sin((2.0* x));
    
    std::cout << "d2/dx2(" << f.NthToString(0)
            << ") = " << f.NthToString(2)
            << " = " << f.Nth(2) <<std::endl;


    return 0;
}
```

**Output**


```
d2/dx2((sin(((CONST[2]) * (VAR[1,ID[2]]))))) = ((CONST[-1]) * (sin(((CONST[2]) * (VAR[1,ID[2]]))))) = -0.909297
```


## Partial Derivatives ##

```

#include <iostream>
#include "ADNumber.hpp"

using namespace std;

typedef ad::ADNumber<double> var;

/*
 * Example of partial derivatives.
 * 
 * reference: http://www.wolframalpha.com/input/?i=d%2Fdx+x%5E2+y%5E4%2C+d%2Fdy+x%5E2+y%5E4&lk=3
 */
int main(int argc, char** argv) {

    
    var f;
    var x = 1;
    var y = 1;
    
    f =  pow(x,2.0)*pow(y,4.0);
    
    std::cout << "{d/dx,d/dy}(" << f.NthToString(0)
            << ") = \n{" << f.NthPartialToString(x,1)<<",\n" << f.NthPartialToString(y,1)<<"}\n"
            << " = {" << f.NthPartial(x,1) <<"," << f.NthPartial(y,1)<<"}" <<std::endl;


    return 0;
}
```

**Output**

```
{d/dx,d/dy}(((pow((VAR[1,ID[2]]),(CONST[2]))) * (pow((VAR[1,ID[3]]),(CONST[4]))))) = 
{(((((CONST[1]) * (CONST[2])) * (pow((VAR[1,ID[2]]),((CONST[2]) - (CONST[1]))))) * (pow((VAR[1,ID[3]]),(CONST[4])))) + ((pow((VAR[1,ID[2]]),(CONST[2]))) * (CONST[0]))),
(((CONST[0]) * (pow((VAR[1,ID[3]]),(CONST[4])))) + ((pow((VAR[1,ID[2]]),(CONST[2]))) * (((CONST[1]) * (CONST[4])) * (pow((VAR[1,ID[3]]),((CONST[4]) - (CONST[1])))))))}
 = {2,4}
```

## Hessian Matrix ##

```
#include <iostream>
#include "ADNumber.hpp"

using namespace std;

typedef ad::ADNumber<double> var;

/*
 * Example of a Hessian matrix.
 * 
 */
int main(int argc, char** argv) {

    double hessian[4];

    var f;
    var x = 1;
    var y = 1;

    f = pow(x, 2.0) * pow(y, 4.0);

    hessian[0] = f.WRT(x, x);
    hessian[1] = f.WRT(x, y);
    hessian[2] = f.WRT(y, x);
    hessian[3] = f.WRT(y, y);

    std::cout << "H(" << f.NthToString(0)
            << ") = \n{[d2f/d2x,d2f/dxdy],[d2f/d2y,d2f/dydx]} =\n"
            <<"{["<<hessian[0]<<","<<hessian[1]<<"],["<<hessian[2]
            <<","<<hessian[3]<<"]}"<<std::endl;


    return 0;
}
```

**Output**

```
H(((pow((VAR[1,ID[2]]),(CONST[2]))) * (pow((VAR[1,ID[3]]),(CONST[4]))))) = 
{[d2f/d2x,d2f/dxdy],[d2f/d2y,d2f/dydx]} =
{[2,8],[8,12]}
```


# ADNumber #

## Public Member Functions ##

```

 ADNumber ()
 ADNumber (T value, T derivative=T(1.0))
 ADNumber (std::string name, T value)
 ADNumber (const ADNumber &orig)
virtual ~ADNumber ()
operator ADNumber< T > () const
operator T() const 
ADNumber< T > & 	operator= (const ADNumber< T > &val)
ADNumber< T > & 	operator= (const T &val)
ADNumber< T > 	operator+ (const ADNumber< T > &rhs) const
ADNumber< T > 	operator+ (const T &rhs) const
ADNumber< T > 	operator- (const ADNumber< T > &rhs) const
ADNumber< T > 	operator- (const T &rhs) const
ADNumber< T > 	operator* (const ADNumber< T > &rhs) const
ADNumber< T > 	operator* (const T &rhs) const
ADNumber< T > 	operator/ (const ADNumber< T > &rhs) const
ADNumber< T > 	operator/ (const T &rhs) const
ADNumber< T > 	operator+= (const ADNumber< T > &rhs)
ADNumber< T > 	operator-= (const ADNumber< T > &rhs)
ADNumber< T > 	operator*= (const ADNumber< T > &rhs)
ADNumber< T > 	operator/= (const ADNumber< T > &rhs)
ADNumber< T > 	operator+= (const T &rhs)
ADNumber< T > 	operator-= (const T &rhs)
ADNumber< T > 	operator*= (const T &rhs)
ADNumber< T > 	operator/= (const T &rhs)
ADNumber< T > 	operator++ ()
ADNumber< T > 	operator-- ()
ADNumber< T > 	operator++ (int)
ADNumber< T > 	operator-- (int)
const T 	GetValue () const
const ADNumber< T > 	WRT (const ADNumber< T > &var0)
const ADNumber< T > 	WRT (const ADNumber< T > &var0, const ADNumber< T > &var1)
const ADNumber< T >	WRT (const ADNumber< T > &var0, const ADNumber< T > &var1, const ADNumber< T > &var2)
const ADNumber< T > 	WRT (const ADNumber< T > &var0, const ADNumber< T > &var1, const ADNumber< T > &var2, const ADNumber< T > &var3)
const ADNumber< T > 	WRT (const ADNumber< T > &var0, const ADNumber< T > &var1, const ADNumber< T > &var2, const ADNumber< T > &var3, const ADNumber< T > &var4)
const ADNumber< T > 	WRT (const ADNumber< T > &var0, const ADNumber< T > &var1, const ADNumber< T > &var2, const ADNumber< T > &var3, const ADNumber< T > &var4, const ADNumber< T > &var5)
const ADNumber< T > 	WRT (const ADNumber< T > &var0, const ADNumber< T > &var1, const ADNumber< T > &var2, const ADNumber< T > &var3, const ADNumber< T > &var4, const ADNumber< T > &var5, const ADNumber< T > &var6)
const ADNumber< T > 	WRT (const ADNumber< T > &var0, const ADNumber< T > &var1, const ADNumber< T > &var2, const ADNumber< T > &var3, const ADNumber< T > &var4, const ADNumber< T > &var5, const ADNumber< T > &var6, const ADNumber< T > &var7)
const ADNumber< T > 	WRT (const ADNumber< T > &var0, const ADNumber< T > &var1, const ADNumber< T > &var2, const ADNumber< T > &var3, const ADNumber< T > &var4, const ADNumber< T > &var5, const ADNumber< T > &var6, const ADNumber< T > &var7, const ADNumber< T > &var8)
const ADNumber< T > 	WRT (const ADNumber< T > &var0, const ADNumber< T > &var1, const ADNumber< T > &var2, const ADNumber< T > &var3, const ADNumber< T > &var4, const ADNumber< T > &var5, const ADNumber< T > &var6, const ADNumber< T > &var7, const ADNumber< T > &var8, const ADNumber< T > &var9)
const ADNumber< T > 	WRT (const std::vector< ADNumber< T > * > &vars)
const std::string 	GetName () const
void 	SetName (const std::string &name)
const uint32_t 	GetID () const
```


## Math Functions ##
```


template<class T >
ad::ADNumber< T > 	atan (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	atan2 (const ad::ADNumber< T > &lhs, const ad::ADNumber< T > &rhs)
template<class T >
ad::ADNumber< T > 	atan2 (T lhs, const ad::ADNumber< T > &rhs)
template<class T >
ad::ADNumber< T > 	atan2 (const ad::ADNumber< T > &lhs, T rhs)
template<class T >
ad::ADNumber< T > 	cos (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	exp (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	log (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	log10 (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	pow (const ad::ADNumber< T > &lhs, const ad::ADNumber< T > &rhs)
template<class T >
ad::ADNumber< T > 	pow (T lhs, const ad::ADNumber< T > &rhs)
template<class T >
ad::ADNumber< T > 	pow (const ad::ADNumber< T > &lhs, T rhs)
template<class T >
ad::ADNumber< T > 	sin (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	sqrt (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	tan (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	acos (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	asin (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	sinh (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	cosh (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	tanh (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	fabs (const ad::ADNumber< T > &val)
template<class T >
ad::ADNumber< T > 	floor (const ad::ADNumber< T > &val)
```