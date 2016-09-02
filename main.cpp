#include "OptimalFilter.h"
/**
* optimal filter design
*
* This project is a part of my private assement of two frameworks:
* 1) the Ultimate++ c++ framework, assesed as a rapid GUI prototyper.
* 2) Eigen linear algebra framework, assesed as a Matlab replacement.
* 
* Anyway, given cutoff frequency (bandwidth; 3db frequency), sampling frequency and filter order,
* OptimalFilter will return an all-pole, near-linear phase low pass filter with optimized
* magnitude response in the passband region.
*
* remarks:
* > since the filter is "near linear" in phase, it should be forward/backward filtered
*   on implementation if non linear phase shift is not acceptable.
* > left click with the mouse on the graphic output to open a label displaying 
*   point information.
*
* Dan I. Malta
**/

/**
* phase unwrapping
**/
void unwrap(Vector<double>& p) {
	int j, N = p.GetCount();
	Vector<double> dp, dps, dp_corr, cumsum;
	
	dp.SetCount(N);
	dps.SetCount(N);
	dp_corr.SetCount(N);
	cumsum.SetCount(N);

    // incremental phase variation 
    for (j = 0; j < N-1; j++) {
        dp[j] = p[j+1] - p[j];
    }
      
    // equivalent phase variation in [-pi, pi]
    for (j = 0; j < N-1; j++) {
        dps[j] = (dp[j] + pi) - floor((dp[j] + pi) / (2 * pi)) * ( 2 * pi) - pi;
    }

    // preserve variation sign for +pi vs. -pi
    for (j = 0; j < N-1; j++) {
        if ((dps[j] == -pi) && (dp[j] > 0)) {
            dps[j] = pi;
        }
    }

    // incremental phase correction
    for (j = 0; j < N-1; j++) {
        dp_corr[j] = dps[j] - dp[j];
    }
      
    // Ignore correction when incremental variation is smaller than cutoff
    for (j = 0; j < N-1; j++) {
        if (fabs(dp[j]) < pi) {
            dp_corr[j] = 0;
        }
    }

    // Find cumulative sum of deltas
    cumsum[0] = dp_corr[0];
    for (j = 1; j < N-1; j++) {
        cumsum[j] = cumsum[j-1] + dp_corr[j];
    }

    // Integrate corrections and add to P to produce smoothed phase values
    for (j = 1; j < N; j++) {
        p[j] += cumsum[j-1];
    }
  }


/**
* optimal low pass FIR design using Durrani-Chapman method
**/
void OptimalFilter::OptimalFilterCalc(const int xi_order,	// filter order
				       const double xi_fc,	// filter cutoff / bandwidth [hz]
				       const double xi_fs) {	// sampling frequency [Hz]
 	// filter denoinator
    // constants
    static const double eps = 1e-9;    
    int k, l, n, j, ind = 0;
    
    // housekeeping
    double fnot = xi_fc / xi_fs;
    int N = xi_order;
    double nn = static_cast<double>(N);
        
    // create the Discrete Prolate Spheroidal Sequences
    MatrixXd sigma(N, N);
    sigma.setZero();
    for (k = 0; k < N - 1; k++) {
        double kk = static_cast<double>(k);
        sigma(k, k) = (nn + 1.0 - 2.0 * (kk + 1.0)) * (nn + 1.0 - 2.0 * (kk + 1.0)) / 4.0 * cos(2 * pi * fnot);
        sigma(k, k + 1) = 0.5 * (kk + 1.0) * (nn - 1.0 - ((kk + 1.0) - 1.0));
        sigma(k + 1, k) = 0.5 * (kk + 1.0) * (nn - (kk + 1.0));
    }
    sigma(N - 1, N - 1) = (nn + 1 - 2 * nn) * (nn + 1 - 2 * nn) / 4 * cos(2 * pi * fnot);
       
    // Discrete Prolate Spheroidal Sequences eigenvalues & eigenvectors
    SelfAdjointEigenSolver<MatrixXd> eigensolver(sigma);
    MatrixXd vv = eigensolver.eigenvectors();	// each column is an eignevector
    ArrayXd lambda = eigensolver.eigenvalues();
    
    // remove zero eigenvalues
    if ((lambda.abs() < eps).any() == 1) {
        for (k = 0; k < lambda.size(); k++) {
            if (abs(lambda(k)) < eps) {
                int kk = lambda.size() - k - 1;
				lambda.segment(k, kk) = lambda.tail(kk);
				lambda.conservativeResize(lambda.size() - 1);
            }
        }
    } 

	// find smallest eigenvalue index
	double minLambda = lambda.minCoeff();
	for (k = 0; k < lambda.size(); k++) {
		if (abs(lambda(k) - minLambda) < 1e-6) {
			ind = k;
			break;
		}
	}

	// extract eigenvecor of smallest eigenvale and normalize it by its first member
	VectorXd v = vv.col(ind);
	v /= v(1);
			
	// filter gain & arguments
	std::complex<double> coeff(0, 2.0 * pi * fnot);		 
	MatrixXcd R(N, N);
	R.setZero();
	for (k = 0; k <= N - 1; k++) {
		for (l = 0; l <= N - 1; l++) {
			double kk = static_cast<double>(k),
				   ll = static_cast<double>(l);
			R(k, l) = std::exp(coeff * ((kk + 1.0) - (ll + 1.0)));
		}
	}
	VectorXcd Arg = v.transpose() * R * v;
	double K = sqrt(1 / MAX2(Arg(0).real(), eps));
	
	// create polynom denominator (convolution between v & itself)
	VectorXd denom(2 * v.size() - 1);
	for (n = 0; n < 2 * v.size() - 1; n++) {
		denom(n) = 0.0;
        for (k = 0; k < v.size(); k++) {
            denom(n) += (k < v.size() ? v[k] : 0.0) *
                        (((n - k < v.size()) && n - k >= 0) ? (v[n - k]) : (0.0));
        }
	}
	denom *= K * K;
	denom /= (denom(0) >= 0) ? (1.0) : (-1.0);
	denom(v.size() - 1) -= 1.0;
	if (abs(denom(0)) > eps) {
		denom /= denom(0);
	} else {
		denom /= eps;
	}

	// find denominator roots (using companion matrix)
	PolynomialSolver<double, Eigen::Dynamic> psolve(denom);
	VectorXcd fact = psolve.roots();
		
	// maintain only stable poles (inside the unit circle)
	factStable.Clear();
	ArrayXcd stablepoles;
	n = 0;
	for (k = 0; k < fact.size(); ++k) {
		if (sqrt(fact(k).real() * fact(k).real() + fact(k).imag() * fact(k).imag()) < 1.0) {
			if ((abs(fact(k).real()) < 1.0) && (abs(fact(k).imag())) < 1.0) {
				factStable.Add(fact(k).real(), fact(k).imag());
				stablepoles.conservativeResize(n+1);
				stablepoles(n) = fact(k);
				n++;
			}
		}
	}

	// locate static gain
	VectorXcd H(1001);
	VectorXcd E(1001);
	std::complex<double> c(0, 2.0 * pi);	
	H.setOnes();
	for (k = 0; k < E.size(); k++) {
		E(k) = std::exp(c * 0.001 * static_cast<double>(k));
	}
	for (k = 0; k < stablepoles.size(); k++) {
		for (j = 0; j < E.size(); j++) {
			H(j) *= E(j) - stablepoles(k);
		}
	}
	
	// transfer function
	for (k = 0; k < H.size(); k++) {
		H(k) = 1.0 / H(k);
	}
	double G = 1.0 / abs(H(0));	// denominator
	H *= G;
	
	// numerator
	VectorXcd polynomial(stablepoles.size() + 1);
  	roots_to_monicPolynomial(stablepoles, polynomial);
  	polynomial.reverseInPlace();
  	
	// filter structutr output
	String s;
	s << "Y[n] = " << G << "*X[n] - ";
	for (k = 1; k < polynomial.size(); k++) {
		s << polynomial(k).real() << "*Y[n-" << k << "] - ";
	}
	s.Remove(s.GetCount() - 2, 2);
	filterOut.SetText(s);
		
	// plot filter poles
	zero.Clear();
	zero.Add(G, 0);
	scatterPoleZero.RemoveAllSeries();
	scatterPoleZero.AddSeries(&circleDraw, 20).NoMark().Stroke(1, Black());//.Dash(LINE_DOTTED);
	scatterPoleZero.AddSeries(zero).MarkStyle<CircleMarkPlot>().Stroke(0, Blue()).MarkColor(Red());
	scatterPoleZero.AddSeries(factStable).MarkStyle<XMarkPlot>().Stroke(0, Red()).MarkColor(Blue());
	xMin = MIN2(MIN2(fact.real().minCoeff(), G) - 0.1, -1.1);
	xMax = MAX2(MAX2(fact.real().maxCoeff(), G) + 0.1, 1.1);
	yMin = MIN2(fact.imag().minCoeff() - 0.1, -1.1);
	yMax = MAX2(fact.imag().maxCoeff() + 0.1, 1.1);
	scatterPoleZero.SetXYMin(xMin, yMin);
	scatterPoleZero.SetRange(abs(xMin) + xMax, abs(yMin) + yMax);
	//scatterPoleZero.SetXYMin(-1.5, -1.5).SetRange(3.0, 3.0);
	
	// calculate filter frequency reponse
	magnitude.Clear();
	phase.Clear();
	Fs.Clear();
	magMin = 999999999.0;
	magMax = -999999999.0;
	//magCross.Clear();
	bool found = false;
	for (k = 0; k < H.size(); ++k) {
		magnitude.Add(20.0 * log10(abs(H(k))));
		phase.Add(r2d * atan2(H(k).imag(), H(k).real()));
		Fs.Add(static_cast<double>(k) * 0.001 * xi_fs);
		
		if (magnitude[k] > magMax) {
			magMax = magnitude[k];
		}
		if (magnitude[k] < magMin) {
			magMin = magnitude[k];
		}
		
		/*
		if ((abs(Fs[k] - xi_fc) <= 0.01) && (found == false)) {
			magCross.Add(xi_fc, magnitude[k]);
			found = true;
		}
		*/
		
		if (Fs[k] >= xi_fs / 2.0) {
			break;
		}
	}
	
	// unwrap phase vector
	phaseMin = 999999999.0;
	phaseMax = -999999999.0;
	unwrap(phase);
	for (k = 0; k < phase.GetCount(); ++k) {
		if (phase[k] > phaseMax) {
			phaseMax = phase[k];
		}
		if (phase[k] < phaseMin) {
			phaseMin = phase[k];
		}
	}
	
	// plot filter frequency responce
	scatterFrequency.RemoveAllSeries();
	scatterFrequency.AddSeries(Fs, magnitude).NoMark().Stroke(1, Blue()).Legend("Magnitude");
	scatterFrequency.AddSeries(Fs, phase).NoMark().Stroke(1, Red()).Legend("Phase").SetDataPrimaryY(false);
	//scatterFrequency.AddSeries(magCross).MarkStyle<CircleMarkPlot>().Stroke(0, Blue()).MarkColor(Blue());
	scatterFrequency.SetDrawY2Reticle();
	scatterFrequency.SetXYMin(0.0, magMin, phaseMin);
	scatterFrequency.SetRange(xi_fs / 2.0,
	                          abs(magMin) + magMax + 0.1,
	                          abs(phaseMin) + phaseMax);

	// calculate step response
	stepMax = -999999999.0;
	stepMin =  999999999.0;
	step.Clear();
	found = false;			// stop calculation flag
	int counter = 0;		// number of while loop calculations
	double u = 0.0,			// filter input (0 @ first cycle, 1 otherwise)
		   ytemp = 0.0,		// filter output (current cycle)
		   dt = 1.0 / xi_fs,	// time step [sec]
	       timeStep = 0.0,			// time index
	       stop = 10.0;			// stoppage criteria
	Vector<double> y;		// filter output vector
	y.Add(0.0);	// initial value
	while (found == false) {
		// filter output
		ytemp = G * u;
		for (k = 1; k < MIN2(polynomial.size(), counter); ++k) {
			if (counter - k > 0) {
				ytemp -= polynomial(k).real() * y[counter - (k - 1)];
			}
		}
		y.Add(ytemp);
		
		// min/max 
		if (ytemp > stepMax) {
			stepMax = ytemp;
		}
		if (ytemp < stepMin) {
			stepMin = ytemp;
		}
		
		// enter step data
		step.Add(timeStep, ytemp);
		
		// stoppage criteria (<2 * filter length> cycles whos sum difference from u is 0.1)
		if (counter > 0) {
			stop = 0.0;
		}
		for (k = 0; k < MIN2(2 * polynomial.size(), counter); ++k) {
			stop += abs(1.0 - y[counter - k]);
		}
		if ((stop <= 0.1) || (counter * dt > 60)) {
			found = true;
		}
		
		// input
		u = 1.0;
		
		// counter's update
		timeStep += dt;
		counter++;
	}
	scatterStepImpulse.AddSeries(step).NoMark().Stroke(1, Blue());
	scatterStepImpulse.SetXYMin(0.0, stepMin - 0.1);
	scatterStepImpulse.SetRange(timeStep, abs(stepMin) + stepMax + 0.1);
}

/**
* calculate
**/
void OptimalFilter::calculateBtnFunc() {
	// retrieve cutoff frequency value & type
	double cutoff = cutoffFrequencyEdit.Scan(cutoffFrequencyEdit.GetText());
	if (cutoffFrequencyType.GetIndex() == 1) {	// transform from [rad/sec] to [hz]
		cutoff *= 0.15915494309189533576888376337251; // = 1 / (2*pi)
	}
	
	// retrieve sampling frequency & type
	double samplingF = samplingFrequencyEdit.Scan(samplingFrequencyEdit.GetText());
	if (samplingFrequencyType.GetIndex() == 1) {	// trgansform from [rad/sec] to [hz]
		samplingF *= 0.15915494309189533576888376337251;	// = 1 / (2*pi)
	}
	
	// do not perform calculation if cutoff frequency is smaller then sampling frequency
	if (samplingF <= cutoff) {
		PromptOK("Cutoff frequency can not be larger or equal to sampling frequency.");
	} else if (samplingF <= 0.0) {
		PromptOK("Sampling frequency can not be negative or zero.");
		samplingFrequencyEdit.SetText("1.0");
	} else {	// every thing is good
		// retrieve filter order
		int orderF = filterOrderSpiner.GetData();
		
		// calculate optimal filter
		OptimalFilterCalc(orderF, cutoff, samplingF);
	}
}

/**
* filter order spinner
**/
void OptimalFilter::filterOrderSpinerFunc() {
	filterOrderEdit.SetText(FormatInt((int)filterOrderSpiner.GetData()));
}

/**
* about dialog
**/
void OptimalFilter::aboutBtnFunc() {
	WithAboutDialogLayout<TopWindow> dlg;
	CtrlLayoutOK(dlg, "About");
	dlg.CenterScreen();
	dlg.Sizeable();
	dlg.Run();
}
	
/**
* constructor
**/
OptimalFilter::OptimalFilter() {
	// load main layout
	CtrlLayout(*this, "Optimal Filter Design Tool");
	
	cutoffFrequencyEdit.SetText(FormatInt((int)10));
	cutoffFrequencyType.Add(0, "Hz")
					   .Add(1, "rad/sec");
	cutoffFrequencyType.SetIndex(0);
	cutoffFrequencyType.Activate();
	
	samplingFrequencyEdit.SetText(FormatInt((int)100));
	samplingFrequencyType.Add(0, "Hz")
					     .Add(1, "rad/sec");
	samplingFrequencyType.SetIndex(0);
	samplingFrequencyType.Activate();
	
	filterOrderSpiner <<= THISBACK(filterOrderSpinerFunc);
	filterOrderSpiner.Step(1);
	filterOrderSpiner.MinMax(2, 20);
	filterOrderSpiner.SetData(2);
	filterOrderEdit.SetText("2");
	
	calculateBtn  <<= THISBACK(calculateBtnFunc);
	aboutBtn      <<= THISBACK(aboutBtnFunc);
	
	// no context menu
	scatterPoleZero.ShowContextMenu(false);
	scatterFrequency.ShowContextMenu(false);
	scatterStepImpulse.ShowContextMenu(false);
	
	// scatter handling
	scatterPoleZero.SetMode(ScatterDraw::MD_ANTIALIASED);
	scatterFrequency.SetMode(ScatterDraw::MD_ANTIALIASED);
	scatterStepImpulse.SetMode(ScatterDraw::MD_ANTIALIASED);
}

/**
* Main
**/
GUI_APP_MAIN {
	OptimalFilter().Run();
}
