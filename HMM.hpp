#ifndef HMM_HPP
#define HMM_HPP

class HMM
{
public:
	// HMM Variables
	int T;								// Sequence Size
	int N;								// Number of states
	vector<double> inits;				// Initial state probabilities
	vector<double> trans;				// Transition probabilities
	vector<double> emits;				// Emission probabilities
	// Viterbi Training Variables
	// Baum Welch Training Variables
	vector<double> alpha;				// Forward probabilities (matrix T*N stored in a vector)
	vector<double> beta;				// Backward probabilities (matrix T*N stored in a vector)
	vector<double> gamma;				// State Occupation Probability (matrix T*N stored in a vector)
	// Storing the transpose of the matrices to fasten the lookup


	// Member Functions
	HMM(int n_states, vector<); 					// Constructor
	~HMM();						 		// Destructor

	// Returns the optimal Viterbi Decode seq in the input referenced vector 
	vector<int>& ViterbiDecode(vector<int>&);
	// Returns the optimal Viterbi Decode seq as a vector of int
	vector<int> ViterbiDecode();

	// Train on the given sequence using Baum Welch algorithm
	void BaumWelchTrain();
	// Train on the given sequence using Viterbi Training algorithm
	void ViterbiTrain();

};

#endif

