#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

/** define user function here:
 *
 * Please following the two steps to write user function:
 * 1. In first step, take x[] and data[] input, output futureState[][] which will be used in interpolation
 * 2. In second step, take x[], data[], vFuture[] and gFuture[] as input, output v, grad[]
 *  */
USER_FUNC_HEAD
{
	/**
	 * preprocessor, do not edit
	 */
	USER_FUNC_PRE;

	/**
	 * step 1:
	 * input: x[], data[]
	 * output: stateFuture[][] of length [nvec, stateDim]
	 */
	// private input ebp
	int POPCOUNT = 0;
	// input stack goes from top to bottom
#define POP(var) double var = data[POPCOUNT++]
#define POPARRAY(var, length) double* var = &data[POPCOUNT]; POPCOUNT+=length

	// individual data
	POP(sigma);
	POP(chi);
	POP(tauL);
	POP(tauR);
	POP(r);
	POP(w);
	POP(tax);

	POP(epsilon);
	POP(a);
	POP(budget);

	// controls
	double ap = x[0];

	double E = budget - ap;
	double u;
	double c;
	double l;
	double inf = 1e20;
	if (E < 1e-12)
		u = -inf;
	else {
		if (E <= (1 - tauL)*w*epsilon / (1 - chi)) {
			l = (1 - chi)*E / ((1 - tauL)*w*epsilon);
			c = chi*E;
		}
		else {
			l = 1;
			c = E - (1 - tauL)*w*epsilon;
		}
	}
	u = pow(c, chi*(1 - sigma))*pow(l, (1 - chi)*(1 - sigma)) / (1 - sigma);

	g_rhs_return[0] = ap;
	g_rhs_return[1] = 1 - l;
	g_rhs_return[2] = c;
	g_rhs_return[3] = u;

	// future State
	stateFuture[0][0] = ap;

	/**
	 * interpolate future state, do not edit
	 */
	USER_FUNC_INTERP;

	/**
	 * step 2:
	 * input: x[], data[], vFuture[], gFuture[][]
	 * output:
	 *    v: scalar value eavaluated at x
	 *    grad[] of size controlDim: gradient of v evaluated at x
	 *    cons[] of size nonlin: evaluation of constraint
	 *    consgrad[][] of size [nonlin, controlDim]: jacobian of cons
	 */
	if (f) {
		v = u + vFuture[0];
	}

	if (grad) {
	}

	if (cons) {
	}

	if (consgradRaw) {
	}

	USER_FUNC_RETURN;
}

#undef MIN
#undef MAX
