#include <string>
#include <fstream>
#include "assert.h"

#ifndef _ANN_H_
#define _ANN_H_

//#define ANN_INLINES
#define ANN_DEBUG_MODE

typedef enum ActivationType
{
  ANN_TANH,
  ANN_LOGISTIC
};

typedef enum RegularizationType
{
  W_DECAY,
  W_ELIMINATION,
  W_APPROXIMATE_SMOOTHER
};

/**
 * The class provides functions for backpropagation learning as well.
 * @brief A class implementing a Multilayer Perceptron with one input, one hidden and one output layer.
 */
class ANN
{
  public:
    ANN(unsigned inputNodes, unsigned hiddenNodes, unsigned outputNodes, unsigned seed = 0, double rndBound = 0.5 /*ANN::RANDOM_BOUND*/);
    ~ANN();

    void initLearnData(double** learnData, double** expectedOutput, unsigned samples, bool makeCopy = false);
    void initValidationData(double** validationData, double** validationExpected, unsigned validationSamples, bool makeCopy = false);
    void learn(unsigned maxIterations, double minError, unsigned seed = 0);
    void learn2(unsigned maxIterations, double minError, unsigned seed = 0);
    void writeParameters(std::string& filename);
    void readParameters(std::string& filename);

    const double* getActivation(double* input);

    std::string getStats() const;

    void useLiteratureValues();
    void normalizeData();
    void normalize(double* sample);
#ifndef ANN_INLINES
    void setFuncType(ActivationType type);
    void setLearningConstant(double value);
    void setMomentumConstant(double value);
    void setValidationMaxIt(unsigned validationMaxIt);
    void setMaxIt(unsigned maxIt);
    void setMinError(double minError);
    void setA(double value);
    void setB(double value);
    void setShuffleSet(bool shuffleSet);
    void setAdjustFrequency(unsigned adjustFrequency);
    void setAdjustConstants(bool adjustConstants);
    void setAdjustFactor(double adjustFactor);
    void setAdjustMin(double adjustMin);
    void setAdjustMax(double adjustMax);
    void setRhoPlus(double rhoPlus);
    void setRhoMinus(double rhoMinus);
    void setLambda(double lambda);
    void setWZero(double w_zero);
    void setRegularization(RegularizationType regularization);
    void setPower(unsigned power);
    void setPeriod(unsigned period);

    unsigned getInputNodes() const;
    unsigned getHiddenNodes() const;
    unsigned getOutputNodes() const;
    ActivationType getFuncType() const;
    double getLearningConstant() const;
    double getMomentumConstant() const;
    unsigned getValidationMaxIt() const;
    unsigned getMaxIt() const;
    double getMinError() const;
    double getA() const;
    double getB() const;
    bool getShuffleSet() const;
    unsigned getAdjustFrequency() const;
    bool getAdjustConstants() const;
    double getAdjustFactor() const;
    double getAdjustMin() const;
    double getAdjustMax() const;
    double getRhoPlus() const;
    double getRhoMinus() const;
    double getLambda() const;
    double getWZero() const;
    RegularizationType getRegularization() const;
    unsigned getPower() const;
    unsigned getPeriod() const;
#else
    inline void setFuncType(ActivationType type) {this->funcType = type;}
    inline void setLearningConstant(double value) {assert(value >= 0);this->learningConstant = value;}
    inline void setMomentumConstant(double value) {assert(value >= 0);this->momentumConstant = value;}
    inline void setValidationMaxIt(unsigned validationMaxIt) {this->validationMaxIt = validationMaxIt;}
    inline void setMaxIt(unsigned maxIt) {this->maxIt = maxIt;}
    inline void setMinError(double minError) {assert(minError >= 0);this->minError = minError;}
    inline void setA(double value) {this->a = value;}
    inline void setB(double value) {this->b = value;}
    inline void setShuffleSet(bool shuffleSet) {this->shuffleSet = shuffleSet;}
    inline void setAdjustFrequency(unsigned adjustFrequency) {this->adjustFrequency = adjustFrequency;}
    inline void setAdjustConstants(bool adjustConstants) {this->adjustConstants = adjustConstants;}
    inline void setAdjustFactor(double adjustFactor) {assert(adjustFactor >= 0);this->adjustFactor = adjustFactor;}
    inline void setAdjustMin(double adjustMin) {assert(adjustMin >= 0);this->adjustMin = adjustMin;}
    inline void setAdjustMax(double adjustMax) {assert(adjustMax >= 0);this->adjustMax = adjustMax;}
    inline void setRhoPlus(double rhoPlus) {assert(rhoPlus > 0);this->rhoPlus = rhoPlus;}
    inline void setRhoMinus(double rhoMinus) {assert(rhoMinus > 0);this->rhoMinus = rhoMinus;}
    inline void setLambda(double lambda) {this->lambda = lambda;}
    inline void setWZero(double w_zero) {this->w_zero = w_zero;}
    inline void setRegularization(RegularizationType regularization) {this->regularization = regularization;}
    inline void setPower(unsigned power) {this->power = power;}
#ifdef ANN_DEBUG_MODE
    inline void setPeriod(unsigned period) {this->period = period;}
#endif

    inline unsigned getInputNodes() const {return this->inputNodes;}
    inline unsigned getHiddenNodes() const {return this->hiddenNodes;}
    inline unsigned getOutputNodes() const {return this->outputNodes;}
    inline ActivationType getFuncType() const {return this->funcType;}
    inline double getLearningConstant() const {return this->learningConstant;}
    inline double getMomentumConstant() const {return this->momentumConstant;}
    inline unsigned getValidationMaxIt() const {return this->validationMaxIt;}
    inline unsigned getMaxIt() const {return this->maxIt;}
    inline double getMinError() const {return this->minError;}
    inline double getA() const {return this->a;}
    inline double getB() const {return this->b;}
    inline bool getShuffleSet() const {return this->shuffleSet;}
    inline unsigned getAdjustFrequency() const {return this->adjustFrequency;}
    inline bool getAdjustConstants() const {return this->adjustConstants;}
    inline double getAdjustFactor() const {return this->adjustFactor;}
    inline double getAdjustMin() const {return this->adjustMin;}
    inline double getAdjustMax() const {return this->adjustMax;}
    inline double getRhoPlus() const {return this->rhoPlus;}
    inline double getRhoMinus() const {return this->rhoMinus;}
    inline double getLambda() const {return this->lambda;}
    inline double getWZero() const {return this->w_zero;}
    inline RegularizationType getRegularization() const {return this->regularization;}
    inline unsigned getPower() const {return this->power;}
#ifdef ANN_DEBUG_MODE
    inline unsigned getPeriod() const {return this->period;}
#endif
#endif

    /** @brief The type of activation function recommended by literature. */
    static const ActivationType LITERATURE_FUNC_TYPE = ANN_TANH;
    /** @brief The recommended maximum amplitude for the recommended activation function. */
    static const double LITERATURE_A_VALUE = 1.7159;
    /** @brief The recommended slope parameter for the recommended activation function. */
    static const double LITERATURE_B_VALUE = 2.0/3.0;
    /** @brief The recommended positive expectation value to use together with the recommended activation function. */
    static const double LITERATURE_POSITIVE_EXPECTED = 1.0;
    /** @brief The recommended negative expectation value to use together with the recommended activation function. */
    static const double LITERATURE_NEGATIVE_EXPECTED = -1.0;
    /** @brief The prefix string for intermediate parameters written out for debugging. */
    static const std::string WRITE_PREFIX_STRING;
    /** @brief The suffix string for intermediate parameters written out for debugging. */
    static const std::string WRITE_SUFFIX_STRING;
  protected:
    double activationFunction(double x);
    double getError(bool validation = false);
    void updateActivation(double* iVec);
    void updateValidation();
    void freeDataset(double** inputData, double** outputData, unsigned setSize);
    std::string memSizeStr(unsigned memSize) const;
    void adjustLearnConstant(unsigned iterations);
    double getPenality() const;
    double getPenalityDerivative(bool output, unsigned i, unsigned j) const;

    /** @brief The number of nodes in the input layer. */
    unsigned inputNodes;
    /** @brief The number of nodes in the hidden layer. */
    unsigned hiddenNodes;
    /** @brief The number of nodes in the output layer. */
    unsigned outputNodes;
    /** @brief The activation of the input layer. */
    double* inputActivation;
    /** @brief The activation of the hidden layer. */
    double* hiddenActivation;
    /** @brief The activation of the output layer. */
    double* outputActivation;
    /** @brief The weight parameters of the hidden layer. */
    double** hiddenParameters;
    /** @brief The weight parameters of the output layer. */
    double** outputParameters;
    /** @brief The type of the activation function of the neurons. */
    ActivationType funcType;
    /** @brief The learning constant used for the learning algorithm. */
    double learningConstant;
    /** @brief The momentum constant used for the learning algorithm. */
    double momentumConstant;
    /** @brief The learning dataset. */
    double** learnData;
    /** @brief The expected outputs for the learning dataset. */
    double** expectedOutput;
    /** @brief The number of samples in the learning dataset. */
    unsigned samples;
    /** @brief Determines whether the learning dataset is a copy or the original learning dataset itself. */
    bool learnDataCopied;

    /** @brief Determines whether validation is done. */
    bool validating;
    /** @brief The number of successive iterations where the current error is above the best error measured so far. */
    unsigned validationIt;
    /** @brief The validation dataset. */
    double** validationData;
    /** @brief The expected values for the validation dataset. */
    double** validationExpected;
    /** @brief The number of samples in the learning dataset. */
    unsigned validationSamples;
    /** @brief Determines whether the validtion dataset is a copy or the original learning dataset itself. */
    bool validationCopied;
    /** @brief The maximum value <code>validationIt</code> can reach before the learning algorithm stops. */
    unsigned validationMaxIt;
    /** @brief The parameters for the hidden neurons that produced the besr error value so far. */
    double** bestHiddenParameters;
    /** @brief The parameters for the output neurons that produced the besr error value so far. */
    double** bestOutputParameters;
    /** @brief The best error value measured so far. */
    double bestErrorValue;

    /** @brief First parameter of the activation function. */
    double a;
    /** @brief Second parameter of the activation function (if needed). */
    double b;
    /** @brief Determines whether the order of the learning dataset should be randomly shuffled between two epochs. */
    bool shuffleSet;

    /** @brief Determines the frequency (in iterations) with which the learning constant should be adjusted. */
    unsigned adjustFrequency;
    /** @brief Determines whether the learning constant should be adjusted while learning. */
    bool adjustConstants;
    /** @brief The factor with which the learning constant should be multiplied when adjusting it during learning. */
    double adjustFactor;
    /** @brief The minimum value the learning constant can not fall below while adjusting it. */
    double adjustMin;
    /** @brief The maximum value the learning constant can not rise above while adjusting it. */
    double adjustMax;
    /** @brief The multiplicative factor to use when increasing the learning constant using <code>ANN::learn2()</code>. */
    double rhoPlus;
    /** @brief The multiplicative factor to use when decreasing the learning constant using <code>ANN::learn2()</code>. */
    double rhoMinus;

    /** @brief Variable representing the importance of the regularization term used. */
    double lambda;
    /** @brief The w0 parameter used for the weight elimination regularization term. */
    double w_zero;
    /** @brief The regularization type used by this <code>ANN</code>. */
    RegularizationType regularization;
    /** @brief The power parameter used for the approximate smoother regularization term. */
    unsigned power;

    /** @brief The default value for the <code>funcType</code> member variable. */
    static const ActivationType DEFAULT_FUNCTION_TYPE = ANN_LOGISTIC;
    /** @brief The default value for the <code>learningConstant</code> member variable. */
    static const double DEFAULT_LEARNING_CONSTANT = 0.1;
    /** @brief The default value for the <code>momentumConstant</code> member variable. */
    static const double DEFAULT_MOMENTUM_CONSTANT = 0;
    /** @brief The bound of the absolute value of the parameters when initializing them randomly. */
    static const double RANDOM_BOUND = 0.5;
    /** @brief The default value for the <code>validationMaxIt</code> member variable. */
    static const unsigned DEFAULT_VALIDATION_MAX_IT = 10;
    /** @brief The default value for the <code>maxIt</code> member variable. */
    static const unsigned DEFAULT_MAX_IT = 1000;
    /** @brief The default value for the <code>minError</code> member variable. */
    static const double DEFAULT_MIN_ERROR = 1e-5;
    /** @brief The default value for the <code>a</code> member variable. */
    static const double DEFAULT_A = 1.0;
    /** @brief The default value for the <code>b</code> member variable. */
    static const double DEFAULT_B = 1.0;
    /** @brief The default value for the <code>shuffleSet</code> member variable. */
    static const bool DEFAULT_SHUFFLE_SET = false;
    /** @brief The default value for the <code>adjustFrequency</code> member variable. */
    static const unsigned DEFAULT_ADJUST_FREQUENCY = 100;
    /** @brief The default value for the <code>adjustConstants</code> member variable. */
    static const bool DEFAULT_ADJUST_CONSTANTS = false;
    /** @brief The default value for the <code>adjustFactor</code> member variable. */
    static const double DEFAULT_ADJUST_FACTOR = 0.1;
    /** @brief The default value for the <code>adjustMin</code> member variable. */
    static const double DEFAULT_ADJUST_MIN = 1e-6;
    /** @brief The default value for the <code>adjustMin</code> member variable. */
    static const double DEFAULT_ADJUST_MAX = 50;
    /** @brief The default value for the <code>rhoMinus</code> member variable. */
    static const double DEFAULT_RHO_MINUS = 0.5;
    /** @brief The default value for the <code>rhoPlus</code> member variable. */
    static const double DEFAULT_RHO_PLUS = 1.2;
    /** @brief The default value for the <code>lambda</code> member variable. */
    static const double DEFAULT_LAMBDA = 1e-10;
    /** @brief The default value for the <code>w_zero</code> member variable. */
    static const double DEFAULT_W_ZERO = 1e-3;
    /** @brief The default value for thy <code>power</code> member variable. */
    static const unsigned DEFAULT_POWER = 2;
    /** @brief The default value for the <code>regularization</code> member variable. */
    static const RegularizationType DEFAULT_REGULARIZATION = W_ELIMINATION;

    /** @brief The maximum number of iterations for the learning algorithm */
    unsigned maxIt;
    /** @brief The minimum error representing the stopping point of the learning algorithm */
    double minError;
#ifdef ANN_DEBUG_MODE
    /** @brief A file for debug output. */
    std::ofstream debugFile;
    /** @brief A file for additional debug output per iteration. (CSV format) */
    std::ofstream dataFile;
    /** @brief The number of iterations after which an intermediate parameter set should be written out. */
    unsigned period;
#endif
};

#endif
