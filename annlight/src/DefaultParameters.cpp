#include "DefaultParameters.h"
#include "SigmoidFunction.h"
#include <iostream>

const std::string DefaultParameters::OPTIMIZATION_GRADIENT_STR = "grad";
const std::string DefaultParameters::OPTIMIZATION_RPROP_STR = "rprop";
const std::string DefaultParameters::OPTIMIZATION_BOTH_STR = "combined";

const std::string DefaultParameters::SF_STEP_STR = "step";
const std::string DefaultParameters::SF_LOGISTIC_STR = "logistic";
const std::string DefaultParameters::SF_LINEAR_STR = "linear";
const std::string DefaultParameters::SF_ARCTAN_STR = "arctan";

const std::string DefaultParameters::STEP_FUNC_ARCTAN_STR = "arctan";
const std::string DefaultParameters::STEP_FUNC_ARCTAN_ZERO_STR = "arctan_zero";
const std::string DefaultParameters::STEP_FUNC_LOGISTIC_STR = "logistic";
const std::string DefaultParameters::STEP_FUNC_TANHYP_STR = "tanh";

const std::string DefaultParameters::FIRST_STAGE_STR = "first";
const std::string DefaultParameters::SECOND_STAGE_STR = "second";

const std::string DefaultParameters::LSM_METHOD_STR = "lsm";
const std::string DefaultParameters::SOFTSTEP_METHOD_STR = "softstep";
const std::string DefaultParameters::HYBRID_METHOD_STR = "hybrid";

const std::string DefaultParameters::GRADIENT_LEARN_PARAM_NAME = "gradient.step";
const std::string DefaultParameters::GRADIENT_LEARN_PARAM_DEFAULT = "0.001";
const std::string DefaultParameters::GRADIENT_PRECISION_NAME = "gradient.precision";
const std::string DefaultParameters::GRADIENT_PRECISION_DEFAULT = "1.0e-3";
const std::string DefaultParameters::RPROP_LEARN_PARAM_NAME = "rprop.step";
const std::string DefaultParameters::RPROP_LEARN_PARAM_DEFAULT = "0.001";
const std::string DefaultParameters::RPROP_PRECISION_NAME = "rprop.precision";
const std::string DefaultParameters::RPROP_PRECISION_DEFAULT = "1.0e-6";
const std::string DefaultParameters::WINDOW_SIZE_NAME = "data.window";
const std::string DefaultParameters::WINDOW_SIZE_DEFAULT = "";
const std::string DefaultParameters::LAMBDA_QUAD_NAME = "lambdaQuad";
const std::string DefaultParameters::LAMBDA_QUAD_DEFAULT = "1e-10";
const std::string DefaultParameters::B_QUAD_NAME = "bQuad";
const std::string DefaultParameters::B_QUAD_DEFAULT = "0";
const std::string DefaultParameters::LAMBDA_LIN_NAME = "lambdaLin";
const std::string DefaultParameters::LAMBDA_LIN_DEFAULT = "1e-10";
const std::string DefaultParameters::B_LIN_NAME = "bLin";
const std::string DefaultParameters::B_LIN_DEFAULT = "0";
const std::string DefaultParameters::LAMBDA_CONST_NAME = "lambdaConst";
const std::string DefaultParameters::LAMBDA_CONST_DEFAULT = "0";
const std::string DefaultParameters::B_CONST_NAME = "bConst";
const std::string DefaultParameters::B_CONST_DEFAULT = "0";
const std::string DefaultParameters::LAMBDA_ALPHA_NAME = "lambdaAlpha";
const std::string DefaultParameters::LAMBDA_ALPHA_DEFAULT = "0.1";
const std::string DefaultParameters::LAMBDA_W1_NAME = "lambdaW1";
const std::string DefaultParameters::LAMBDA_W1_DEFAULT = "0.1";
const std::string DefaultParameters::LAMBDA_W2_NAME = "lambdaW2";
const std::string DefaultParameters::LAMBDA_W2_DEFAULT = "0.1";
const std::string DefaultParameters::LAMBDA_BETA_NAME = "lambdaBeta";
const std::string DefaultParameters::LAMBDA_BETA_DEFAULT = "0.1";
const std::string DefaultParameters::OPT_MODE_NAME = "optMode";
const std::string DefaultParameters::OPT_MODE_DEFAULT = DefaultParameters::OPTIMIZATION_GRADIENT_STR;
const std::string DefaultParameters::ALPHA_PLUS_NAME = "alpha_plus";
const std::string DefaultParameters::ALPHA_PLUS_DEFAULT = "";
const std::string DefaultParameters::BETA_PLUS_NAME = "beta_plus";
const std::string DefaultParameters::BETA_PLUS_DEFAULT = "";
const std::string DefaultParameters::GAMMA_PLUS_NAME = "gamma_plus";
const std::string DefaultParameters::GAMMA_PLUS_DEFAULT = "";
const std::string DefaultParameters::ALPHA_MINUS_NAME = "alpha_minus";
const std::string DefaultParameters::ALPHA_MINUS_DEFAULT = "";
const std::string DefaultParameters::BETA_MINUS_NAME = "beta_minus";
const std::string DefaultParameters::BETA_MINUS_DEFAULT = "";
const std::string DefaultParameters::GAMMA_MINUS_NAME = "gamma_minus";
const std::string DefaultParameters::GAMMA_MINUS_DEFAULT = "";
const std::string DefaultParameters::W_PLUS_NAME = "w_minus";
const std::string DefaultParameters::W_PLUS_DEFAULT = "0.5";
const std::string DefaultParameters::W_MINUS_NAME = "w_plus";
const std::string DefaultParameters::W_MINUS_DEFAULT = "0.5";
const std::string DefaultParameters::QUADRATIC_NAME = "quadratic";
const std::string DefaultParameters::QUADRATIC_DEFAULT = "0";
const std::string DefaultParameters::LINEAR_NAME = "linear";
const std::string DefaultParameters::LINEAR_DEFAULT = "1";
const std::string DefaultParameters::OPTIMIZER_SEED_NAME = "optimizer.seed";
const std::string DefaultParameters::OPTIMIZER_SEED_DEFAULT = "0";
const std::string DefaultParameters::OPTIMIZER_RANDOM_INIT_NAME = "optimizer.rinit";
const std::string DefaultParameters::OPTIMIZER_RANDOM_INIT_DEFAULT = "1";
const std::string DefaultParameters::OPTIMIZER_SFUNC_TYPE_NAME = "optimizer.func_type";
const std::string DefaultParameters::OPTIMIZER_SFUNC_TYPE_DEFAULT = "arctan";
const std::string DefaultParameters::LOG_ALPHA_NAME = "sf.alpha";
const std::string DefaultParameters::LOG_ALPHA_DEFAULT = "1";
const std::string DefaultParameters::LOG_TYPE_NAME = "sf.type";
const std::string DefaultParameters::LOG_TYPE_DEFAULT = DefaultParameters::SF_LOGISTIC_STR;
const std::string DefaultParameters::RPROP_RHOPLUS_NAME = "rprop.rhoplus";
const std::string DefaultParameters::RPROP_RHOPLUS_DEFAULT = "1.5";
const std::string DefaultParameters::RPROP_RHOMINUS_NAME = "rprop.rhominus";
const std::string DefaultParameters::RPROP_RHOMINUS_DEFAULT = "0.5";
const std::string DefaultParameters::RPROP_WMIN_NAME = "rprop.wmin";
const std::string DefaultParameters::RPROP_WMIN_DEFAULT = "0.000001";
const std::string DefaultParameters::RPROP_WMAX_NAME = "rprop.wmax";
const std::string DefaultParameters::RPROP_WMAX_DEFAULT = "10";
const std::string DefaultParameters::GRADIENT_MAX_ITERATIONS_NAME = "gradient.maxit";
const std::string DefaultParameters::GRADIENT_MAX_ITERATIONS_DEFAULT = "500";
const std::string DefaultParameters::RPROP_MAX_ITERATIONS_NAME = "rprop.maxit";
const std::string DefaultParameters::RPROP_MAX_ITERATIONS_DEFAULT = "500";
const std::string DefaultParameters::DEBUGLOGGING_NAME = "logging.debug";
const std::string DefaultParameters::DEBUGLOGGING_DEFAULT = "0";
const std::string DefaultParameters::DATALOGGING_NAME = "logging.data";
const std::string DefaultParameters::DATALOGGING_DEFAULT = "0";
const std::string DefaultParameters::LEARN_METHOD_NAME = "learn_method";
const std::string DefaultParameters::LEARN_METHOD_DEFAULT = DefaultParameters::SOFTSTEP_METHOD_STR;
const std::string DefaultParameters::VALIDATION_MAXIT_NAME = "validation.maxit";
const std::string DefaultParameters::VALIDATION_MAXIT_DEFAULT = "10";
const std::string DefaultParameters::OPTIMIZER_STAGE_NAME = "optimizer.stage";
const std::string DefaultParameters::OPTIMIZER_STAGE_DEFAULT = DefaultParameters::FIRST_STAGE_STR;
const std::string DefaultParameters::OPTIMIZE_W_NAME = "optimizeW";
const std::string DefaultParameters::OPTIMIZE_W_DEFAULT = "0";
const std::string DefaultParameters::WPLUSMAX_NAME = "w_plus_max";
const std::string DefaultParameters::WPLUSMAX_DEFAULT = "0.95";
const std::string DefaultParameters::WPLUSMIN_NAME = "w_plus_min";
const std::string DefaultParameters::WPLUSMIN_DEFAULT = "0.05";
const std::string DefaultParameters::WREGPOWER_NAME = "wRegPower";
const std::string DefaultParameters::WREGPOWER_DEFAULT = "1";
const std::string DefaultParameters::OPTIMIZE_ALPHA_NAME = "optimizeAlpha";
const std::string DefaultParameters::OPTIMIZE_ALPHA_DEFAULT = "0";
const std::string DefaultParameters::ALPHAMAX_NAME = "alpha_max";
const std::string DefaultParameters::ALPHAMAX_DEFAULT = "10";
const std::string DefaultParameters::ALPHAMIN_NAME = "alpha_min";
const std::string DefaultParameters::ALPHAMIN_DEFAULT = "0";
const std::string DefaultParameters::OPTIMIZE_BETA_NAME = "optimizeBeta";
const std::string DefaultParameters::OPTIMIZE_BETA_DEFAULT = "0";
const std::string DefaultParameters::BETAMAX_NAME = "beta_max";
const std::string DefaultParameters::BETAMAX_DEFAULT = "10";
const std::string DefaultParameters::BETAMIN_NAME = "beta_min";
const std::string DefaultParameters::BETAMIN_DEFAULT = "0";
const std::string DefaultParameters::INPUT_ENTRIES_NAME = "inputEntries";
const std::string DefaultParameters::INPUT_ENTRIES_DEFAULT = "0";
const std::string DefaultParameters::CYCLIC_RESTART_WRITE_NAME = "cyclicRestartWrite" ;
const std::string DefaultParameters::CYCLIC_RESTART_WRITE_DEFAULT = "0";

/**
 * @brief Standard constructor for the <code>Parameters</code> class.
 */
DefaultParameters::DefaultParameters()
{
  this->parameterList.push_back(make_pair(DefaultParameters::GRADIENT_LEARN_PARAM_NAME, DefaultParameters::GRADIENT_LEARN_PARAM_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::GRADIENT_PRECISION_NAME, DefaultParameters::GRADIENT_PRECISION_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::RPROP_LEARN_PARAM_NAME, DefaultParameters::RPROP_LEARN_PARAM_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::RPROP_PRECISION_NAME, DefaultParameters::RPROP_PRECISION_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::WINDOW_SIZE_NAME, DefaultParameters::WINDOW_SIZE_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LAMBDA_QUAD_NAME, DefaultParameters::LAMBDA_QUAD_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::B_QUAD_NAME, DefaultParameters::B_QUAD_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LAMBDA_LIN_NAME, DefaultParameters::LAMBDA_LIN_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::B_LIN_NAME, DefaultParameters::B_LIN_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LAMBDA_CONST_NAME, DefaultParameters::LAMBDA_CONST_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::B_CONST_NAME, DefaultParameters::B_CONST_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LAMBDA_ALPHA_NAME, DefaultParameters::LAMBDA_ALPHA_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LAMBDA_W1_NAME, DefaultParameters::LAMBDA_W2_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LAMBDA_W2_NAME, DefaultParameters::LAMBDA_W2_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LAMBDA_BETA_NAME, DefaultParameters::LAMBDA_BETA_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::OPT_MODE_NAME, DefaultParameters::OPT_MODE_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::ALPHA_PLUS_NAME, DefaultParameters::ALPHA_PLUS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::BETA_PLUS_NAME, DefaultParameters::BETA_PLUS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::GAMMA_PLUS_NAME, DefaultParameters::GAMMA_PLUS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::ALPHA_MINUS_NAME, DefaultParameters::ALPHA_MINUS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::BETA_MINUS_NAME, DefaultParameters::BETA_MINUS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::GAMMA_MINUS_NAME, DefaultParameters::GAMMA_MINUS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::W_PLUS_NAME, DefaultParameters::W_PLUS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::W_MINUS_NAME, DefaultParameters::W_MINUS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::QUADRATIC_NAME, DefaultParameters::QUADRATIC_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LINEAR_NAME, DefaultParameters::LINEAR_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::OPTIMIZER_SEED_NAME, DefaultParameters::OPTIMIZER_SEED_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::OPTIMIZER_RANDOM_INIT_NAME, DefaultParameters::OPTIMIZER_RANDOM_INIT_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LOG_ALPHA_NAME, DefaultParameters::LOG_ALPHA_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LOG_TYPE_NAME, DefaultParameters::LOG_TYPE_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::RPROP_RHOPLUS_NAME, DefaultParameters::RPROP_RHOPLUS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::RPROP_RHOMINUS_NAME, DefaultParameters::RPROP_RHOMINUS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::RPROP_WMIN_NAME, DefaultParameters::RPROP_WMIN_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::RPROP_WMAX_NAME, DefaultParameters::RPROP_WMAX_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::RPROP_MAX_ITERATIONS_NAME, DefaultParameters::RPROP_MAX_ITERATIONS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::GRADIENT_MAX_ITERATIONS_NAME, DefaultParameters::GRADIENT_MAX_ITERATIONS_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::DEBUGLOGGING_NAME, DefaultParameters::DEBUGLOGGING_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::DATALOGGING_NAME, DefaultParameters::DATALOGGING_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::LEARN_METHOD_NAME, DefaultParameters::LEARN_METHOD_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::OPTIMIZER_SFUNC_TYPE_NAME, DefaultParameters::OPTIMIZER_SFUNC_TYPE_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::VALIDATION_MAXIT_NAME, DefaultParameters::VALIDATION_MAXIT_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::OPTIMIZER_STAGE_NAME, DefaultParameters::OPTIMIZER_STAGE_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::OPTIMIZE_W_NAME, DefaultParameters::OPTIMIZE_W_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::WPLUSMAX_NAME, DefaultParameters::WPLUSMAX_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::WPLUSMIN_NAME, DefaultParameters::WPLUSMIN_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::WREGPOWER_NAME, DefaultParameters::WREGPOWER_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::OPTIMIZE_ALPHA_NAME, DefaultParameters::OPTIMIZE_ALPHA_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::ALPHAMAX_NAME, DefaultParameters::ALPHAMAX_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::ALPHAMIN_NAME, DefaultParameters::ALPHAMIN_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::OPTIMIZE_BETA_NAME, DefaultParameters::OPTIMIZE_BETA_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::BETAMAX_NAME, DefaultParameters::BETAMAX_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::BETAMIN_NAME, DefaultParameters::BETAMIN_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::INPUT_ENTRIES_NAME, DefaultParameters::INPUT_ENTRIES_DEFAULT));
  this->parameterList.push_back(make_pair(DefaultParameters::CYCLIC_RESTART_WRITE_NAME, DefaultParameters::CYCLIC_RESTART_WRITE_DEFAULT));
}

/**
 * @brief Standard destructor for the <code>Parameters</code> class.
 */
DefaultParameters::~DefaultParameters()
{
}

/**
 * This variable contains pairs of strings where the first element is the name of a parameter,
 * and the second element is the parameter's default value.
 * @brief Returns the <code>parameterList</code> member variable.
 *
 * @return The <code>parameterList</code> member variable
 */
const std::list<std::pair<std::string, std::string> >& DefaultParameters::getParameters()
{
  return this->parameterList;
}

/**
 * @brief Returns a string representation of a <code>OptimizationScheme</code> variable.
 *
 * @param scheme A <code>OptimizationScheme</code> that has to be transformed into a string
 * @return A string representation of the given <code>OptimizationScheme</code>
 */
std::string DefaultParameters::toString(OptimizationScheme scheme)
{
  switch (scheme)
  {
    case GRADIENT:
      return DefaultParameters::OPTIMIZATION_GRADIENT_STR;
    case RPROP:
      return DefaultParameters::OPTIMIZATION_RPROP_STR;
    case COMBINED:
      return DefaultParameters::OPTIMIZATION_BOTH_STR;
    default:
      std::cerr << "Unknown optimization scheme: " << scheme << std::endl;
      exit(1);
  }
}

/**
 * @brief Returns a string representation of a Sigmoid function type.
 *
 * @param type The Sigmoid function type as string
 * @return A string representation of the given Sigmoid function type
 */
std::string DefaultParameters::typeToString(short type)
{
  switch (type)
  {
    case SigmoidFunction::LOGISTIC:
      return DefaultParameters::SF_LOGISTIC_STR;
    case SigmoidFunction::LINEAR:
      return DefaultParameters::SF_LINEAR_STR;
    case SigmoidFunction::ARCTAN:
      return DefaultParameters::SF_ARCTAN_STR;
    case SigmoidFunction::STEP:
      return DefaultParameters::SF_STEP_STR;
    default:
      std::cerr << "Unknown step function type: " << type << std::endl;
      exit(1);
  }
}

/**
 * @brief Returns a short value representing a Sigmoid function type.
 *
 * @param str The string representation of a Sigmoid function type
 */
short DefaultParameters::toSFType(std::string str)
{
  if (str == DefaultParameters::SF_LOGISTIC_STR)
    return SigmoidFunction::LOGISTIC;
  else if (str == DefaultParameters::SF_LINEAR_STR)
    return SigmoidFunction::LINEAR;
  else if (str == DefaultParameters::SF_STEP_STR)
    return SigmoidFunction::STEP;
  else if (str == DefaultParameters::SF_ARCTAN_STR)
    return SigmoidFunction::ARCTAN;
  else
  {
    std::cerr << "Unknown step function type string: " << str << std::endl;
    exit(1);
  }
}

/**
 * @brief Returns the <code>SoftStepFuncType</code> represented by the given string parameter.
 *
 * @param str The string representation of a <code>SoftStepFuncType</code>
 * @return The <code>SoftStepFuncType</code> represented by <code>str</code>
 */
SoftStepFuncType DefaultParameters::toSoftStepFuncType(std::string str)
{
  if (str == DefaultParameters::STEP_FUNC_ARCTAN_STR)
    return ARCTAN;
  else if (str == DefaultParameters::STEP_FUNC_ARCTAN_ZERO_STR)
    return ARCTAN_ZERO;
  else if (str == DefaultParameters::STEP_FUNC_LOGISTIC_STR)
    return LOGISTIC;
  else if (str == DefaultParameters::STEP_FUNC_TANHYP_STR)
    return TANHYP;
  else
  {
    std::cerr << "Unknown soft step function type string:" << str << std::endl;
    exit(1);
  }
}

/**
 * @brief Returns the <code>OptimizationScheme</code> represented by the given string parameter.
 *
 * @param str The string representation of a <code>OptimizationScheme</code>
 * @return The <code>OptimizationScheme</code> represented by <code>str</code>
 */
OptimizationScheme DefaultParameters::toOptimizationScheme(std::string str)
{
  if (str == DefaultParameters::OPTIMIZATION_GRADIENT_STR)
    return GRADIENT;
  else if (str == DefaultParameters::OPTIMIZATION_RPROP_STR)
    return RPROP;
  else if (str == DefaultParameters::OPTIMIZATION_BOTH_STR)
    return COMBINED;
  else
  {
    std::cerr << "Unknown optimization scheme string: " << str << std::endl;
    exit(1);
  }
}

/**
 * @brief Returns the <code>StageType</code> represented by the given string parameter.
 *
 * @param str The string representation of a <code>StageType</code>
 * @return The <code>StageType</code> represented by <code>str</code>
 */
StageType DefaultParameters::toStageType(std::string str)
{
  if (str == FIRST_STAGE_STR)
    return FIRST;
  else if (str == SECOND_STAGE_STR)
    return SECOND;
  else
  {
    std::cerr << "Unknown stage type string: " << str << std::endl;
    exit(1);
  }
}

/**
 * @brief Returns the <code>LearnMethodType</code> represented by the given string parameter.
 *
 * @param str The string representation of a <code>LearnMethodType</code>
 * @return The <code>LearnMethodType</code> representated by <code>str</code>
 */
LearnMethodType DefaultParameters::toLearnMethodType(std::string str)
{
  if (str == LSM_METHOD_STR)
    return LSM;
  else if (str == SOFTSTEP_METHOD_STR)
    return SOFTSTEP;
  else if (str == HYBRID_METHOD_STR)
    return HYBRID;
  else
  {
    std::cerr << "Unknown learn method type string: " << str << std::endl;
    exit(1);
  }
}
