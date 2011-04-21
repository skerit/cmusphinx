%feature("docstring", "
    Log-space math class.
    
    This class provides fast logarithmic math functions in base
    1.000+epsilon, useful for fixed point speech recognition.

    @param base: The base B in which computation is to be done.
    @type base: float
    @param shift: Log values are shifted right by this many bits.
    @type shift: int
    @param use_table Whether to use an add table or not
    @type use_table: bool") LogMath;
typedef struct logmath_s {
} LogMath;

%extend LogMath {
	LogMath(double base, int shift=0, bool use_table=false) {
		return logmath_init(base, shift, use_table);
	}
	~LogMath() {
		logmath_free($self);
	}
	%feature("docstring","
        Get the log base.

        @return: Logarithmic base used in this object.
        @rtype: float
") get_base;
	double get_base() {
		return logmath_get_base($self);
	}
	%feature("docstring","
        Get the log-zero value.

        @return: Smallest number representable by this object.
        @rtype: int
") get_zero;
	int get_zero() {
		return logmath_get_zero($self);
	}
	%feature("docstring","
        Add two numbers in log-space without using an add-table.

        @param a: Logarithm A.
        @type a: int
        @param b: Logarithm B.
        @type b: int
        @return: log(exp(a)+exp(b))
        @rtype: int
") add;
	int add_exact(int p, int q) {
		return logmath_add_exact($self, p, q);
	}
	%feature("docstring","
        Add two numbers in log-space.

        @param a: Logarithm A.
        @type a: int
        @param b: Logarithm B.
        @type b: int
        @return: log(exp(a)+exp(b))
        @rtype: int
") add;
	int add(int p, int q) {
		return logmath_add($self, p, q);
	}
	%feature("docstring","
        Return log-value of a number.

        @param x: Number (in linear space)
        @type x: float
        @return: Log-value of x.
        @rtype: int
") log;
	int log(double p) {
		return logmath_log($self, p);
	}
	%feature("docstring", "
        Return linear of a log-value

        @param x: Logarithm X (in this object's base)
        @type x: int
        @return: Exponent (linear value) of X.
        @rtype: float
") exp;
	double exp(int logb_p) {
		return logmath_exp($self, logb_p);
	}
	%feature("docstring", "
        Return log-value of a natural logarithm.

        @param x: Logarithm X (in base e)
        @type x: float
        @return: Log-value equivalent of x.
        @rtype: int
") ln_to_log;
	int ln_to_log(double log_p) {
		return logmath_ln_to_log($self, log_p);
	}
	%feature("docstring", "
        Return natural logarithm of a log-value.

        @param x: Logarithm X (in this object's base)
        @type x: int
        @return: Natural log equivalent of x.
        @rtype: float
") log_to_ln;
	double log_to_ln(int logb_p) {
		return logmath_log_to_ln($self, logb_p);
	}
	%feature("docstring","
        Return log-value of a base 10 logarithm.

        @param x: Logarithm X (in base 10)
        @type x: float
        @return: Log-value equivalent of x.
        @rtype: int
") log10_to_log;
	int log10_to_log(double log_p) {
		return logmath_log10_to_log($self, log_p);
	}
	%feature("docstring","
        Return logarithm in base 10 of a log-value.

        @param x: Logarithm X (in this object's base)
        @type x: int
        @return: log10 equivalent of x.
        @rtype: float
") log_to_log10;
	double log_to_log10(int logb_p) {
		return logmath_log_to_log10($self, logb_p);
	}
};
