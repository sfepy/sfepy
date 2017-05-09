//contres.h
#ifndef CONTRES_H
#define CONTRES_H

#ifdef __cplusplus
extern "C" {  // only need to export C interface if used by C++ source code
#endif

#ifdef __linux__

	void sfd2(double* H, double* dH, double r);

	void sfd6(double* H, double* dH, double r, double s);

	void getAABB(double* AABBmin, double* AABBmax, int nsd, int nnod, double* X, double longestEdge, int* IEN, int* ISN, int* elementID, int* segmentID, int n, int nsn, int nes, int nen, int neq);

  void assembleContactResidualAndStiffness(double* Gc, double* vals, double* rows, double* cols, int* len, double* GPs, int* ISN, int* IEN, double* X, double* U, double* H, double* dH, double* gw, double* activeGPsOld, int neq, int nsd, int npd, int ngp, int nes, int nsn, int nen, int GPs_len, double epss, bool keyContactDetection, bool keyAssembleKc);

	void getLongestEdgeAndGPs(double* longestEdge, double* GPs, int n, int nsd, int ngp, int neq, int nsn, int nes, int nen, int* elementID, int* segmentID, int* ISN, int* IEN, double* H, double* X);

	void evaluateContactConstraints(double* GPs, int* ISN, int* IEN, int* N, double* AABBmin, double* AABBmax, int* head, int* next, double* X, int* elementID, int* segmentID, int n, int nsn, int nsd, int npd, int ngp, int nen, int nes, int neq, double longestEdge);
#elif _WIN32
			  // windows code goes here

	void __declspec(dllexport) sfd2(double* H, double* dH, double r);
	void __declspec(dllexport) sfd6(double* H, double* dH, double r, double s);
	void __declspec(dllexport) getAABB(double* AABBmin, double* AABBmax, int nsd, int nnod, double* X, double longestEdge, int* IEN, int* ISN, int* elementID, int* segmentID, int n, int nsn, int nes, int nen, int neq);

	void __declspec(dllexport) assembleContactResidualAndStiffness(double* Gc, double* vals, double* rows, double* cols, int* len, double* GPs, int* ISN, int* IEN, double* X, double* U, double* H, double* dH, double* gw, double* activeGPsOld, int neq, int nsd, int npd, int ngp, int nes, int nsn, int nen, int GPs_len, double epss, bool keyContactDetection, bool keyAssembleKc);

	void __declspec(dllexport) getLongestEdgeAndGPs(double* longestEdge, double* GPs, int n, int nsd, int ngp, int neq, int nsn, int nes, int nen, int* elementID, int* segmentID, int* ISN, int* IEN, double* H, double* X);

	void __declspec(dllexport) evaluateContactConstraints(double* GPs, int* ISN, int* IEN, int* N, double* AABBmin, double* AABBmax, int* head, int* next, double* X, int* elementID, int* segmentID, int n, int nsn, int nsd, int npd, int ngp, int nen, int nes, int neq, double longestEdge);
#else

#endif

#ifdef __cplusplus
}
#endif

#endif  // CONTRES_H
