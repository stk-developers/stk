
namespace STK
{
  
  /** **************************************************************************
   * @brief Abstract class for binary question
   */
  class BinaryQuestion
  {
  public:
    typedef bool Answer;
    typedef xxx  Term;

  public:

    /** 
     * @brief Evaluates the term rTerm and gives binary answer
     * @param rTerm Termo to evaluate
     * @return bool
     */
    virtual Answer
    Eval(const Term& rTerm) const = 0;

    virtual BinaryQuestion*
    Clone();

  }; // class BinaryQuestion


  /** **************************************************************************
   * @brief Simple predictor question 
   *
   * The principle of this question is: "is the previous n'th predictor variable
   * such and such???
   */
  class PredictorQuestion : public BinaryQuestion
  {
  private:
    std::vector<bool> mSet;
    int               mPredictor;

  public:
    /** 
     * @brief Simple constructor
     */
    PredictorQuestion();

    /** 
     * @brief Copy constructor
     * @param rOrig 
     */
    PredictorQuestion(const PredictorQuestion& rOrig);

    /** 
     * @brief Virtual destructor
     */
    virtual
    ~PredictorQuestion();


    virtual Answer
    Eval(const Term& rTerm) const;
  }; // class PredictorQuestion

} // namespace STK

//################################# EOF ########################################
