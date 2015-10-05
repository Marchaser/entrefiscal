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
	POP(theta);
	POP(gamma);
	POP(lambda);
	POP(delta);
	POP(tauPi);
	POP(tauR);
	POP(z);
	POP(r);
	POP(w);
	POP(tax);

	POP(zeta);
	POP(a);
	POP(budget);

	// controls
	double ap = x[0];

	// production
	/*
	double A = z*zeta;
	double coef1 = (1 - gamma)*(r + delta) / (gamma*w);
	double k = pow(A*theta*gamma*pow(coef1, theta*(1 - gamma)) / (r + delta), 1 / (1 - theta));
	k = MIN(k, a*lambda);
	double n = coef1*k;
	double pi = A*pow(k, gamma*theta)*pow(n, (1 - gamma)*theta) - (r + delta)*k - w*n;
	*/

	double c = budget - ap;
	double u;
	double cLow = 1e-6;
	double uAtCLow = pow(cLow, chi*(1 - sigma)) / (1 - sigma);
	double upAtCLow = chi*(1 - sigma)*uAtCLow / cLow;
	if (c > cLow)
		u = pow(c, chi*(1 - sigma)) / (1 - sigma);
	else
		u = uAtCLow + upAtCLow*(c - cLow);

	g_rhs_return[0] = ap;
	// g_rhs_return[1] = k;
	// g_rhs_return[2] = n;
	g_rhs_return[1] = c;
	g_rhs_return[2] = u;

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
