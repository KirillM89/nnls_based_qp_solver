#include "gtest/gtest.h"
#include "test_utils.h"
#include "test_data.h"
#include "decorator.h"
using namespace QP_NNLS;
//using namespace QP_NNLS_TEST_DATA;

TEST(Utils, MatrixMult1) {
	matrix_t M1 = { {1.0, 2.0}, {3.0, 4.0} };
	matrix_t baseline = { {7.0, 10.0}, {15.0, 22.0} };
	TestMatrixMult(M1, M1, baseline);
}
TEST(Utils, MatrixMult2) {
	//square 3x3
	matrix_t M1 = { {1.0, 2.0, -1.0}, {3.0, 4.0, -7.0}, {1.0, -5.0, 3.0} };
	matrix_t M2 = { {5.0, 2.0, -7.0}, {1.0, -2.0, -9.0}, {0, 6.0, -4.0} };
	matrix_t baseline = { {7.0, -8.0, -21.0}, {19.0, -44.0, -29.0}, {0.0 ,30.0, 26.0} };
	TestMatrixMult(M1, M2, baseline);
}
TEST(Utils, MatrixMult3) {
	//1x2  2x2
	matrix_t M1 = {{1.0, 2.0}};
	matrix_t M2 = {{3.0, 4.0}, {5.0, 6.0}};
	matrix_t baseline = {{13.0, 16.0}};
	TestMatrixMult(M1, M2, baseline);
}
TEST(Utils, MatrixMult4) {
	// 2x1  1x2
	matrix_t M1 = {{1.0}, {-5.0}};
	matrix_t M2 = {{3.0, -4.0}};
	matrix_t baseline = {{3.0, -4.0}, {-15.0, 20.0}};
	TestMatrixMult(M1, M2, baseline);
}
TEST(Utils, MultTransp1) {
	matrix_t M = {{1.0}, {2.0}};
	std::vector<double> v = {1.0, 2.0};
	std::vector<double> baseline = {5.0};
	TestMatrixMultTranspose(M, v, baseline);
}
TEST(Utils, MultTransp2) {
	matrix_t M = { {-5.0} };
	std::vector<double> v = {5.0};
	TestMatrixMultTranspose(M, v, {-25.0});
}
TEST(Utils, MultTransp3) {
	matrix_t M = { {1.0, 2.0}, {3.0, 4.0} };
	std::vector<double> v = { 5.0, -1.0};
	std::vector<double> baseline = { 2.0, 6.0 };
	TestMatrixMultTranspose(M, v, baseline);
}
TEST(Utils, MultTransp4) {
	matrix_t M = { {1.0, 2.0, 3.0}, {3.0, 4.0, 5.0} };
	std::vector<double> v = { 5.0, -1.0};
	std::vector<double> baseline = { 2.0, 6.0, 10.0 };
	TestMatrixMultTranspose(M, v, baseline);
}
TEST(Utils, MultTransp5) {
	matrix_t M = {{1.0, 2.0, 3.0, -2.0}, {3.0, 4.0, 5.0, 0.0}};
	std::vector<double> v = {5.0, -1.0};
	std::vector<double> baseline = {2.0, 6.0, 10.0, -10.0};
	TestMatrixMultTranspose(M, v, baseline);
}

TEST(Utils, M1M2T_1) {
	const matrix_t M = { {1.0, 2.0}, {3.0, 4.0} };
	const matrix_t baseline = { {5.0, 11.0}, {11.0, 25.0} };
	TestM1M2T(M, M, baseline);
}
TEST(Utils, M1M2T_2) {
	const matrix_t M = { {-1.0, 2.0, -5.0}, {3.0, 4.0, -10.0}, {2.0 , 1.0, -1.0 } };
	const matrix_t baseline = { {30.0, 55.0, 5.0}, {55.0, 125.0, 20.0} ,{5.0, 20.0, 6.0} };
	TestM1M2T(M, M, baseline);
}
TEST(Utils, M1M2T_3) {
	const matrix_t M = {{5.0}, {7.0}};
	const matrix_t baseline = {{25.0, 35.0}, {35.0, 49.0}};
	TestM1M2T(M, M, baseline);
}
TEST(Utils, M1M2T_Simple4) {
	const matrix_t M1 = {{5.0}, {7.0}};
	const matrix_t M2 = {{5.0}, {7.0}, {9.0}};
	const matrix_t baseline = { {25.0, 35.0, 45.0}, {35.0, 49.0, 63.0} };
	TestM1M2T(M1, M2, baseline);
}
TEST(Utils, M1M2T_Simple5) {
	const matrix_t M1 = {{5.0}, {7.0}, {-1.0}};
	const matrix_t M2 = {{2.0}, {-5.0}, {9.0}};
	const matrix_t baseline = {{10.0, -25.0, 45.0}, {14.0, -35.0, 63.0}, {-2.0, 5.0, -9.0}};
	TestM1M2T(M1, M2, baseline);
}
TEST(Utils, M1M2T_Simple6) {
	const matrix_t M1 = {{5.0}};
	const matrix_t M2 = {{1.0}, {5.0}, {9.0}};
	const matrix_t baseline = {{5.0, 25.0, 45.0}};
	TestM1M2T(M1, M2, baseline);
}
TEST(Utils, M1M2T_Simple7) {
	const matrix_t M1 = {{5.0}, {-1.0}, {4.0}};
	const matrix_t M2 = {{1.0}};
	const matrix_t baseline = M1;
	TestM1M2T(M1, M2, baseline);
}
TEST(Utils, M1M2T_Simple8) {
	const matrix_t M1 = {{5.0, 3.0}, {-1.0, 2.0}, {4.0, -1.0}};
	const matrix_t M2 = {{1.0, 2.0}};
	const matrix_t baseline = {{11.0}, {3.0},{2.0}};
	TestM1M2T(M1, M2, baseline);
}

TEST(Utils, InPlaceLDLT_T1) {
    matrix_t M = {{1.0, 0.0, 0.0},
                  {0.0, 1.0, 0.0},
                  {0.0, 0.0, 1.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_T2) {
    matrix_t M = {{2.0, 0.0, 0.0},
                  {0.0, 2.0, 0.0},
                  {0.0, 0.0, 2.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_T3) {
    matrix_t M = {{2.0, 1.0, 0.0},
                  {1.0, 2.0, 0.0},
                  {0.0, 0.0, 2.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_T4) {
    matrix_t M = {{15.0, 10.0, 0.0},
                  {10.0, 2.0, 7.0},
                  {0.0, 7.0, 2.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_T5) {
    matrix_t M = {{1.0, -1.0, 3.0},
                  {-1.0, 1.1, 0.0},
                  {3.0, 0.0, 1.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_T6) {
    matrix_t M = {{1.0, -1.0, 3.0},
                  {-1.0, 1.0001, 0.0},
                  {3.0, 0.0, 1.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_T10) {
    matrix_t M = {{2.0, -1.0, 3.0, 0.0},
                  {-1.0, 2.0, 0.0, 7.0},
                  {3.0, 0.0, -4.0, -5.0},
                  {0.0, 7.0, -5.0, -5.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_T11) {
    matrix_t M = {{1.0, -1.0, 3.0},
                  {-1.0, 2.0, 0.1},
                  {3.0, 0.1, 1.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_Swap2Elmts) {
    matrix_t M = {{1.0 , -1.0, 3.0, 7.0},
                  {-1.0, 1.0 , 0.1, -4.0},
                  {3.0 , 0.1 , 1.0, 1.0},
                  {7.0 ,-4.0 , 1.0, 10.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_1SwapT1) {
    //swap 0 and 3
    matrix_t M = {{1.0 , -1.0, 3.0, 7.0},
                  {-1.0, 1.0 , 0.1, -4.0},
                  {3.0 , 0.1 , 1.0, 1.0},
                  {7.0 ,-4.0 , 1.0, 10.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_1SwapT2) {
    //swap 0,3;
    matrix_t M = {{1.0 , -1.0, 3.0, 7.0},
                  {-1.0, 2.0 , 0.1, -4.0},
                  {3.0 , 0.1 , 1.0, 1.0},
                  {7.0 ,-4.0 , 1.0, 10.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_1SwapT3) {
    //swap 1,2;
    matrix_t M = {{ 5.0 , 0.1, -4.0},
                  { 0.1 , 1.0, 1.0},
                  {-4.0 , 1.0, 3.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_2SwapT1) {
    //swap 0,3; 1,2
    matrix_t M = {{1.0 , -1.0, 3.0, 7.0},
                  {-1.0, 1.0 , 0.1, -4.0},
                  {3.0 , 0.1 , 3.0, 1.0},
                  {7.0 ,-4.0 , 1.0, 10.0}};
    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_PosSmdf1) {
    //swap 0,3; 1,2
    matrix_t M = {{1.0 , 2.0 , 3.0},
                  {2.0 , 2.0 , -2.0},
                  {3.0 ,-2.0 , 6.0}};

    TestInPlaceLdlt(M);
}
TEST(Utils, InPlaceLDLT_HS118) {
    matrix_t M;
    GetIdentityMatrix(15, M);
    for(std::size_t i = 0; i < 15; ++i) {
        if (i == 2 || i == 5 || i == 8 || i == 11 || i == 14) {
            M[i][i] = 0.0003;
        } else {
            M[i][i] = 0.0002;
        }
    }
    TestInPlaceLdlt(M);
}
TEST(Utils_LDLT, Test1) {
    const matrix_t M = { {1.0}, {1.0} };
    const std::vector<double> S = {1.0, 1.0};
    const std::set<unsigned int> active = {0, 1};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Test2) {
    const matrix_t M = { {1.0}, {2.0} };
    const std::vector<double> S = {3.0, 4.0};
    const std::set<unsigned int> active = {0, 1};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Test3) {
    const matrix_t M = {{1.0}};
    const std::vector<double> S = {3.0};
    const std::set<unsigned int> active = {0};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Test4) {
    const matrix_t M = {{1.0, 2.0, 3.0}, {-10.5, 100.7, 23.5}};
    const std::vector<double> S = {3.0, -400.987};
    const std::set<unsigned int> active = {0, 1};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Test5) {
    const matrix_t M = {{1.0, 2.0, 3.0}, {-10.5, 100.7, 23.5}, {0.0, 0.0, 0.0}};
    const std::vector<double> S = {3.0, -400.987, 1.0};
    const std::set<unsigned int> active = {0, 1, 2};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Test6) {
    const matrix_t M = {{1.0, 2.0, 3.0}, {-10.0, -10.0, -10.0}, {0.0, 0.0, 0.0}};
    const std::vector<double> S = {3.0, -10.0, 0.0};
    const std::set<unsigned int> active = {0, 1, 2};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Test7) {
    const matrix_t M = {{1.0, 2.0, 3.0}, {-10.0, -10.0, -10.0}, {5.0, 5.0, 5.0}};
    const std::vector<double> S = {3.0, -10.0, 5.0};
    const std::set<unsigned int> active = {0, 1, 2};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Test8) {
    const matrix_t M = {{1.0, 2.0, 3.0}, {-10.0, -10.0, -10.0}, {5.0, 5.0, 5.0}};
    const std::vector<double> S = {3.0, -10.0, 5.0};
    const std::set<unsigned int> active = {0, 2};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Test9) {
    const matrix_t M = {{1.0, 2.0, 3.0},
                        {-10.0, -10.0, -10.0},
                        {5.0, 5.0, 5.0},
                        {0.002, 10.33, -25.9},
                        {-6.0, 10.0, -500.5 }};
    const std::vector<double> S = {3.0, -10.0, 5.0, 0.0, -9.0};
    const std::set<unsigned int> active = {0, 2, 4};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Test10) {
    const matrix_t M = {{1.0, 2.0, -3.0},
                        {-10.0, -10.0, 10.0},
                        {5.0, 5.0, 5.0},
                        {-0.002, 10.33, 0.9},
                        {-6.0, 10.0, 1500.5 }};
    const std::vector<double> S = {3.0, -10.0, 15.0, 0.0, -9.0};
    const std::set<unsigned int> active = {0, 4};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Test11) {
    const matrix_t M = {{1.0, 2.0, -3.0},
                        {-10.0, -10.0, 10.0},
                        {5.0, 5.0, 5.0},
                        {-0.002, 10.33, 0.9},
                        {-6.0, 10.0, 1500.5 }};
    const std::vector<double> S = {3.0, -10.0, 15.0, 0.0, -9.0};
    std::set<unsigned int> active = {3};
    LDLT ldlt(M, S);
    TestLDLT(ldlt, M, S, active);
    active = {1, 4};
    TestLDLT(ldlt, M, S, active);
    active = {0, 2};
    TestLDLT(ldlt, M, S, active);
    active = {2, 3};
    TestLDLT(ldlt, M, S, active);
    active = {0};
    TestLDLT(ldlt, M, S, active);
}
TEST(Utils_LDLT, Add) {
    const matrix_t M = {{1.0, 2.0, -3.0},
                        {-10.0, -10.0, 10.0},
                        {5.0, 5.0, 5.0},
                        {-0.002, 10.33, 0.9},
                        {-6.0, 10.0, 1500.5 }};
    const std::vector<double> S = {3.0, -10.0, 15.0, 0.0, -9.0};
    LdltTester tester;
    tester.Set(M, S);
    tester.Add(1);
    tester.Add(3);
    tester.Add(2);
    tester.Add(0);
    tester.Add(4);
}
TEST(Utils_LDLT, AddPvt) {
    const matrix_t M = {{1.0, 2.0, -3.0},
                        {-10.0, -10.0, 10.0},
                        {5.0, 5.0, 5.0},
                        {-0.002, 10.33, 0.9},
                        {-6.0, 10.0, 1500.5 }};
    const std::vector<double> S = {3.0, -10.0, 15.0, 0.0, -9.0};
    LdltTester tester;
    tester.Set(M, S);
    tester.AddPvt(1);
    tester.AddPvt(3);
    tester.AddPvt(2);
    tester.AddPvt(0);
    tester.AddPvt(4);
}
TEST(Utils_LDLT, DelPvt1) {
    const matrix_t M = {{1.0, 2.0, -3.0},
                        {-10.0, -10.0, 10.0},
                        {5.0, 5.0, 5.0},
                        {-0.002, 10.33, 0.9},
                        {-6.0, 10.0, 1500.5 }};
    const std::vector<double> S = {3.0, -10.0, 15.0, 0.0, -9.0};
    LdltTester tester;
    tester.Set(M, S);
    //case 1
    tester.AddPvt(1);
    tester.Delete(1);
    //case 2
    tester.AddPvt(2);
    tester.AddPvt(0);
    tester.Delete(0);
    tester.Delete(2);
    //case 3             // matrix size:
    tester.AddPvt(1);    // 1
    tester.AddPvt(4);    // 2
    tester.Delete(1);    // 1
    tester.AddPvt(1);    // 2
    tester.Delete(4);    // 1
    tester.AddPvt(3);    // 2
    tester.AddPvt(0);    // 3
    tester.Delete(3);    // 2
    tester.AddPvt(2);    // 3
}
TEST(Utils_LDLT, DelPvt2) {
    //case with equal norms
    const matrix_t M = {{1.0, 2.0, -3.0},
                        {-1.0, -3.0, 2.0},
                        {5.0, 5.0, 5.0},
                        {-0.002, 10.33, 0.9},
                        {-1.5, 0.0 , 11.5 }};
    const std::vector<double> S = {3.0, -3.0, 15.0, 0.0, 4.0};
    //norms2[0] = 23 norms2[1] = 23
    LdltTester tester;
    tester.Set(M, S);
    //case 1
    tester.AddPvt(1);
    tester.Delete(1);
    //case 2
    tester.AddPvt(2);
    tester.AddPvt(0);
    tester.Delete(0);
    tester.Delete(2);
    //case 3             // matrix size:
    tester.AddPvt(1);    // 1
    tester.AddPvt(4);    // 2
    tester.Delete(1);    // 1
    tester.AddPvt(1);    // 2
    tester.Delete(4);    // 1
    tester.AddPvt(3);    // 2
    tester.AddPvt(0);    // 3
    tester.Delete(3);    // 2
    tester.AddPvt(2);    // 3
}
TEST(Utils_LDLT, Del)
{
    const matrix_t M = {{1.0, 2.0, -3.0},
                        {-10.0, -10.0, 10.0},
                        {5.0, 5.0, 5.0},
                        {-0.002, 10.33, 0.9},
                        {-6.0, 10.0, 1500.5 }};
    const std::vector<double> S = {3.0, -10.0, 15.0, 0.0, -9.0};
    LdltTester tester;
    tester.Set(M, S);
    //case 1
    tester.Add(1);
    tester.Delete(1);
    //case 2
    tester.Add(2);
    tester.Add(0);
    tester.Delete(0);
    tester.Delete(2);
    //case 3          // matrix size:
    tester.Add(1);    // 1
    tester.Add(4);    // 2
    tester.Delete(1); // 1
    tester.Add(1);    // 2
    tester.Delete(4); // 1
    tester.Add(3);    // 2
    tester.Add(0);    // 3
    tester.Delete(3); // 2
    tester.Add(2);    // 3
}
using namespace QP_NNLS_TEST_DATA;
TEST(QP_PROBLEM, SIMPLE_1_NO_BOUNDS)
{
    QPProblem(SIMPLE_1_NO_BOUNDS::in, SIMPLE_1_NO_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_2_NO_BOUNDS)
{
    QPProblem(SIMPLE_2_NO_BOUNDS::in, SIMPLE_2_NO_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_3_NO_BOUNDS)
{
    QPProblem(SIMPLE_3_NO_BOUNDS::in, SIMPLE_3_NO_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_4_NO_BOUNDS)
{
    QPProblem(SIMPLE_4_NO_BOUNDS::in, SIMPLE_4_NO_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_5_NO_BOUNDS)
{
    QPProblem(SIMPLE_5_NO_BOUNDS::in, SIMPLE_5_NO_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_6_NO_BOUNDS)
{
    QPProblem(SIMPLE_6_NO_BOUNDS::in, SIMPLE_6_NO_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_7_NO_BOUNDS)
{
    QPProblem(SIMPLE_7_NO_BOUNDS::in, SIMPLE_7_NO_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_8_NO_BOUNDS)
{
    QPProblem(SIMPLE_8_NO_BOUNDS::in, SIMPLE_8_NO_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_9_NO_BOUNDS)
{
    QPProblem(SIMPLE_9_NO_BOUNDS::in, SIMPLE_9_NO_BOUNDS::out);
}
TEST(DISABLED_QP_PROBLEM, SIMPLE_1_LW_BOUNDS)
{
    QPProblem(SIMPLE_1_LW_BOUNDS::in, SIMPLE_1_LW_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_2_LW_BOUNDS)
{
    QPProblem(SIMPLE_2_LW_BOUNDS::in, SIMPLE_2_LW_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_3_LW_BOUNDS)
{
    QPProblem(SIMPLE_3_LW_BOUNDS::in, SIMPLE_3_LW_BOUNDS::out);
}
TEST(QP_PROBLEM, SIMPLE_1_UP_BOUNDS)
{
    QPProblem(SIMPLE_1_UP_BOUNDS::in, SIMPLE_1_UP_BOUNDS::out);
}








