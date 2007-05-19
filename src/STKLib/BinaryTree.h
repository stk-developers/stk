
namespace STK
{
  class Distribution
  {
  private:
    double                  mN;       ///< soft counts
    std::vector<ProbType>   mDistr;

  public:
    Distribution() : mN(0.0), mDistr() 
    {}

    double
    Entropy() const;

    const size_t
    Size() const;
  };

  /** 
   * @brief Binary tree representation
   */
  class DBTree 
  {
  private:
    DBTree*           mpTree0;
    DBTree*           mpTree1;
    BinaryQuestion*   mpQuestion;

  public:

    const bool
    IsLeaf() const;
  }; // class BinaryTree
}
//##################################### EOF ###################################/
