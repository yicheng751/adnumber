# ADNumber #

## Public Member Functions ##

```

 ADNumber ()
 ADNumber (T value, T derivative=T(1.0))
 ADNumber (std::string name, T value, T derivative=T(1.0))
 ADNumber (const ADNumber &orig)
virtual ~ADNumber ()
operator ADNumber< T > () const
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