#include <list>
#include <string>

#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

/**
 * Enumeration for switching the optimization mode.
 */
typedef enum OptimizationScheme
{
  GRADIENT,  // Use gradient descent procedure only
  RPROP,     // Use RPROP procedure only
  COMBINED   // Use both gradient descent procedure and RPROP procedure
};

typedef enum SoftStepFuncType
{
  ARCTAN,
  ARCTAN_ZERO,
  LOGISTIC,
  TANHYP
};

typedef enum StageType
{
  FIRST,
  SECOND
};

typedef enum LearnMethodType
{
  LSM,
  SOFTSTEP,
  HYBRID
};

/**
 * @brief Class holding all parameters with name and default value.
 */
class DefaultParameters
{
  public:
    DefaultParameters();
    ~DefaultParameters();

    const std::list<std::pair<std::string, std::string> >& getParameters();

    /** @brief The name of the learning step parameter of the gradient based learning procedure. */
    static const std::string GRADIENT_LEARN_PARAM_NAME;
    /** @brief The default value of the learning step parameter of the gradient based learning procedure. */
    static const std::string GRADIENT_LEARN_PARAM_DEFAULT;
    /** @brief The name of the precision parameter of the gradient based learning procedure. */
    static const std::string GRADIENT_PRECISION_NAME;
    /** @brief The default value of the precision parameter of the gradient based learning procedure. */
    static const std::string GRADIENT_PRECISION_DEFAULT;
    /** @brief The name of the learning step parameter of the RPROP learning procedure. */
    static const std::string RPROP_LEARN_PARAM_NAME;
    /** @brief The default value of the learning step parameter of the RPROP learning procedure. */
    static const std::string RPROP_LEARN_PARAM_DEFAULT;
    /** @brief The name of the precision parameter of the RPROP learning algorithm. */
    static const std::string RPROP_PRECISION_NAME;
    /** @brief The default value of the precision parameter of the RPROP learning algorithm. */
    static const std::string RPROP_PRECISION_DEFAULT;
    /** @brief The name of the window size parameter. */
    static const std::string WINDOW_SIZE_NAME;
    /** @brief The default value of the window size parameter. */
    static const std::string WINDOW_SIZE_DEFAULT;
    /** @brief The name of the vector dimension parameter. */
    static const std::string VEC_DIM_NAME;
    /** @brief The default value of the vector dimension parameter. */
    static const std::string VEC_DIM_DEFAULT;
    /** @brief The name of the objective function's lambda parameter for quadratic terms. */
    static const std::string LAMBDA_QUAD_NAME;
    /** @brief The default value of the objective function's lambda parameter for quadratic terms. */
    static const std::string LAMBDA_QUAD_DEFAULT;
    /** @brief The name of the objective function's b parameter for quadratic terms. */
    static const std::string B_QUAD_NAME;
    /** @brief The default value of the objective function's b parameter for quadratic terms. */
    static const std::string B_QUAD_DEFAULT;
    /** @brief The name of the objective function's lambda parameter for linear terms. */
    static const std::string LAMBDA_LIN_NAME;
    /** @brief The default value of the objective function's lambda parameter for linear terms. */
    static const std::string LAMBDA_LIN_DEFAULT;
    /** @brief The name of the objective function's b parameter for linear terms. */
    static const std::string B_LIN_NAME;
    /** @brief The default value of the objective function's b parameter for linear terms. */
    static const std::string B_LIN_DEFAULT;
    /** @brief The name of the objective function's lambda parameter for the constant term. */
    static const std::string LAMBDA_CONST_NAME;
    /** @brief The default value of the objective function's lambda parameter for the constant term. */
    static const std::string LAMBDA_CONST_DEFAULT;
    /** @brief The name of the objective function's b parameter for the constant term. */
    static const std::string B_CONST_DEFAULT;
    /** @brief The default value of the objective function's b parameter for the constant term. */
    static const std::string B_CONST_NAME;
    /** @brief The name of the objective function's lambda parameter for the alpha parameter. */
    static const std::string LAMBDA_ALPHA_NAME;
    /** @brief The default value of the objective function's lambda parameter for the alpha parameter. */
    static const std::string LAMBDA_ALPHA_DEFAULT;
    /** @brief The name of the objective function's lambda parameter for the alpha parameter. */
    static const std::string LAMBDA_BETA_NAME;
    /** @brief The default value of the objective function's lambda parameter for the alpha parameter. */
    static const std::string LAMBDA_BETA_DEFAULT;
    /** @brief The name of the objective function's first lambda parameter for the w parameter. */
    static const std::string LAMBDA_W1_NAME;
    /** @brief The default value of the objective function's lambda parameter for the alpha parameter. */
    static const std::string LAMBDA_W1_DEFAULT;
    /** @brief The name of the objective function's lambda parameter for the w parameter. */
    static const std::string LAMBDA_W2_NAME;
    /** @brief The default value of the objective function's lambda parameter for the alpha parameter. */
    static const std::string LAMBDA_W2_DEFAULT;
    /** @brief The name of the <code>Optimizer</code> optimization scheme parameter. */
    static const std::string OPT_MODE_NAME;
    /** @brief The default value of the <code>Optimizer</code> optimization scheme parameter. */
    static const std::string OPT_MODE_DEFAULT;
    /** @brief The name of the positive step function's alpha parameter. */
    static const std::string ALPHA_PLUS_NAME;
    /** @brief The default value of the positive step function's alpha parameter. */
    static const std::string ALPHA_PLUS_DEFAULT;
    /** @brief The name of the positive step function's beta parameter. */
    static const std::string BETA_PLUS_NAME;
    /** @brief The default value of the positive step function' beta parameter. */
    static const std::string BETA_PLUS_DEFAULT;
    /** @brief The name of the positive step function's gamma parameter. */
    static const std::string GAMMA_PLUS_NAME;
    /** @brief The default value of the positive step function' gamma parameter. */
    static const std::string GAMMA_PLUS_DEFAULT;
    /** @brief The name of the negative step function's alpha parameter. */
    static const std::string ALPHA_MINUS_NAME;
    /** @brief The default value of the negative step function's alpha parameter. */
    static const std::string ALPHA_MINUS_DEFAULT;
    /** @brief The name of the negative step function's beta parameter. */
    static const std::string BETA_MINUS_NAME;
    /** @brief The default value of the negative step function's beta parameter. */
    static const std::string BETA_MINUS_DEFAULT;
    /** @brief The name of the negative step function's gamma parameter. */
    static const std::string GAMMA_MINUS_NAME;
    /** @brief The default value of the negative step function's gamma parameter. */
    static const std::string GAMMA_MINUS_DEFAULT;
    /** @brief The name of the objective function's positive weight parameter. */
    static const std::string W_PLUS_NAME;
    /** @brief The default value of the objective function's positive weight parameter. */
    static const std::string W_PLUS_DEFAULT;
    /** @brief The name of the objective function's negative weight parameter. */
    static const std::string W_MINUS_NAME;
    /** @brief The default value of the objective function's negative weight parameter. */
    static const std::string W_MINUS_DEFAULT;
    /** @brief The name of the <code>Optimizer</code> quadratic term parameter. */
    static const std::string QUADRATIC_NAME;
    /** @brief The default value of the <code>Optimizer</code> quadratic term parameter. */
    static const std::string QUADRATIC_DEFAULT;
    /** @brief The name of the <code>Optimizer</code> linear term parameter */
    static const std::string LINEAR_NAME;
    /** @brief The default value of the <code>Optimizer</code> linear term parameter. */
    static const std::string LINEAR_DEFAULT;
    /** @brief The name of the <code>Optimizer</code> seed parameter. */
    static const std::string OPTIMIZER_SEED_NAME;
    /** @brief The default value of the <code>Optimizer</code> seed parameter. */
    static const std::string OPTIMIZER_SEED_DEFAULT;
    /** @brief The name of the <code>Optimizer</code> random init parameter. */
    static const std::string OPTIMIZER_RANDOM_INIT_NAME;
    /** @brief The default value of the <code>Optimizer</code> random init parameter. */
    static const std::string OPTIMIZER_RANDOM_INIT_DEFAULT;
    /** @brief The name of the <code>SigmoidFunction</code> alpha parameter. */
    static const std::string LOG_ALPHA_NAME;
    /** @brief The default value of the <code>SigmoidFunction</code> alpha parameter. */
    static const std::string LOG_ALPHA_DEFAULT;
    /** @brief The name of the <code>SigmoidFunction</code> type parameter. */
    static const std::string LOG_TYPE_NAME;
    /** @brief The default value of the <code>SigmoidFunction</code> type parameter. */
    static const std::string LOG_TYPE_DEFAULT;
    /** @brief The name of the positive rho parameter of the RPROP algorithm. */
    static const std::string RPROP_RHOPLUS_NAME;
    /** @brief The default value of the positive rho parameter of the RPROP algorithm. */
    static const std::string RPROP_RHOPLUS_DEFAULT;
    /** @brief The name of the negative rho parameter of the RPROP algorithm. */
    static const std::string RPROP_RHOMINUS_NAME;
    /** @brief The default value of the negative rho parameter of the RPROP algorithm. */
    static const std::string RPROP_RHOMINUS_DEFAULT;
    /** @brief The name of the minimum weight parameter of the RPROP algorithm. */
    static const std::string RPROP_WMIN_NAME;
    /** @brief The default value of the minimum weight parameter of the RPROP algorithm. */
    static const std::string RPROP_WMIN_DEFAULT;
    /** @brief The name of the maximum weight parameter of the RPROP algorithm. */
    static const std::string RPROP_WMAX_NAME;
    /** @brief The default value of the maximum weight parameter of the RPROP algorithm. */
    static const std::string RPROP_WMAX_DEFAULT;
    /** @brief The name of the maximum number of iterations parameter for the gradient based learning procedure. */
    static const std::string GRADIENT_MAX_ITERATIONS_NAME;
    /** @brief The default value of the maximum number of iterations parameter for the gradient based learning procedure. */
    static const std::string GRADIENT_MAX_ITERATIONS_DEFAULT;
    /** @brief The name of the maximum number of iterations parameter for the RPROP learning procedure. */
    static const std::string RPROP_MAX_ITERATIONS_NAME;
    /** @brief The default value of the maximum number of iterations parameter for the RPROP learning procedure. */
    static const std::string RPROP_MAX_ITERATIONS_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s debug logging parameter. */
    static const std::string DEBUGLOGGING_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s debug logging parameter. */
    static const std::string DEBUGLOGGING_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s data logging parameter. */
    static const std::string DATALOGGING_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s data logging parameter. */
    static const std::string DATALOGGING_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s use LSM parameter. */
    static const std::string LEARN_METHOD_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s use LSM parameter. */
    static const std::string LEARN_METHOD_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s maximum validation iterations parameter. */
    static const std::string VALIDATION_MAXIT_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s maximum validation iterations parameter. */
    static const std::string VALIDATION_MAXIT_DEFAULT;
    /** @brief The name of the <code>Optimizers</code>'s soft step function type parameter. */
    static const std::string OPTIMIZER_SFUNC_TYPE_NAME;
    /** @brief The default value of the <code>Optimizers</code>'s soft step function type parameter. */
    static const std::string OPTIMIZER_SFUNC_TYPE_DEFAULT;
    /** @brief The name of the <code>Optimizers</code>'s stage type parameter. */
    static const std::string OPTIMIZER_STAGE_NAME;
    /** @brief The default of the <code>Optimizers</code>'s stage type parameter. */
    static const std::string OPTIMIZER_STAGE_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s optimize w parameter. */
    static const std::string OPTIMIZE_W_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s optimize w parameter. */
    static const std::string OPTIMIZE_W_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s <code>w_plus_max</code> parameter. */
    static const std::string WPLUSMAX_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s <code>w_plus_max</code> parameter. */
    static const std::string WPLUSMAX_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s <code>w_plus_min</code> parameter. */
    static const std::string WPLUSMIN_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s <code>w_plus_min</code> parameter. */
    static const std::string WPLUSMIN_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s <code>wRegPower</code> parameter. */
    static const std::string WREGPOWER_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s <code>wRegPower</code> parameter. */
    static const std::string WREGPOWER_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s optimize alpha parameter. */
    static const std::string OPTIMIZE_ALPHA_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s optimize alpha parameter. */
    static const std::string OPTIMIZE_ALPHA_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s maximum alpha parameter value. */
    static const std::string ALPHAMAX_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s maximum alpha parameter value. */
    static const std::string ALPHAMAX_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s minimum alpha parameter value. */
    static const std::string ALPHAMIN_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s alpha parameter value. */
    static const std::string ALPHAMIN_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s optimize beta parameter. */
    static const std::string OPTIMIZE_BETA_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s optimize beta parameter. */
    static const std::string OPTIMIZE_BETA_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s maximum beta parameter value. */
    static const std::string BETAMAX_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s maximum beta parameter value. */
    static const std::string BETAMAX_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s minimum beta parameter value. */
    static const std::string BETAMIN_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s beta parameter value. */
    static const std::string BETAMIN_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s input entries value. */
    static const std::string INPUT_ENTRIES_NAME;
    /** @brief The default value of the <code>Optimizer</code>'s input entries value. */
    static const std::string INPUT_ENTRIES_DEFAULT;
    /** @brief The name of the <code>Optimizer</code>'s cyclic restart writing parameter. */
    static const std::string CYCLIC_RESTART_WRITE_NAME;
    /** @brief The dafault value of the <code>Optimizer</code>'s cyclic restart writing parameter. */
    static const std::string CYCLIC_RESTART_WRITE_DEFAULT;

    static std::string toString(OptimizationScheme scheme);
    static OptimizationScheme toOptimizationScheme(std::string str);
    static SoftStepFuncType toSoftStepFuncType(std::string str);
    static StageType toStageType(std::string str);
    static LearnMethodType toLearnMethodType(std::string str);

    static std::string typeToString(short type);
    static short toSFType(std::string str);

  protected:
    /** @brief A list of all existing parameters with name and default value. */
    std::list<std::pair<std::string, std::string> > parameterList;

    /** @brief String representation of the optimization procedure using gradient descent only. */
    static const std::string OPTIMIZATION_GRADIENT_STR;
    /** @brief String representation of the optimization procedure using RPROP only. */
    static const std::string OPTIMIZATION_RPROP_STR;
    /** @brief String representation of the optimization procedure using both gradient descent and RPROP. */
    static const std::string OPTIMIZATION_BOTH_STR;
    /** @brief String representation of the step Sigmoid function type. */
    static const std::string SF_STEP_STR;
    /** @brief String representation of the logistic Sigmoid function type. */
    static const std::string SF_LOGISTIC_STR;
    /** @brief String representation of the linear Sigmoid function type. */
    static const std::string SF_LINEAR_STR;
    /** @brief String representation of the arcus tangens Sigmoid function type. */
    static const std::string SF_ARCTAN_STR;
    /** @brief String representation of the arcus tangens soft step function type. */
    static const std::string STEP_FUNC_ARCTAN_STR;
    /** @brief String representation of the arcus tangens soft step function type bounded by 0 and 1. */
    static const std::string STEP_FUNC_ARCTAN_ZERO_STR;
    /** @brief String representation of the logistic soft step function type. */
    static const std::string STEP_FUNC_LOGISTIC_STR;
    /** @brief String representation of the tangens hyperbolicus soft step function type. */
    static const std::string STEP_FUNC_TANHYP_STR;
    /** @brief String representation of the first stage type. */
    static const std::string FIRST_STAGE_STR;
    /** @brief String representation of the second stage type. */
    static const std::string SECOND_STAGE_STR;
    /** @brief String representation of the LSM learning method. */
    static const std::string LSM_METHOD_STR;
    /** @brief String representation of the Soft Step function learning method. */
    static const std::string SOFTSTEP_METHOD_STR;
    /** @brief String representation of the hybrid learning method. */
    static const std::string HYBRID_METHOD_STR;
};

#endif
