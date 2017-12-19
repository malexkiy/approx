#include <approx.hpp>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm>
#include <sstream>
#include <iterator>
#include <catch.hpp>

SCENARIO("approx test", "[test]") {
	std::vector<double>
		X = { 0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9 },
		Y = { 3, 4.5, 1.7, 0.7, -1, -2, 4 },
		coefs,
		apY,
		_coefs = { 2.9999999988, 93.6609630735, -1253.5720197526, 5798.561576767, -12500.0828384738, 12581.8454197684, -4747.8506020892 },
		_apY = { 2.9999999998453211, 4.5000000010598793, 1.6999999978024625, 0.70000000162359388, -1.0000000004023093, -1.9999999999242846, 3.9999999999945430 };

	coefs = nma::polyfit(X, Y, 6);
	
	std::ostringstream oss1, oss2;
	std::copy(_coefs.begin(), _coefs.end() - 1,
		std::ostream_iterator<double>(oss1, ","));

	std::copy(coefs.begin(), coefs.end() - 1,
		std::ostream_iterator<double>(oss2, ","));

	oss1 << _coefs.back();
	oss2 << coefs.back();

	std::string s1 = oss1.str(),
		s2 = oss2.str();

	REQUIRE(s1 == s2);

	
	
	apY = nma::polyval(coefs, X);

	std::ostringstream oss3, oss4;
	std::copy(_apY.begin(), _apY.end() - 1,
		std::ostream_iterator<double>(oss3, ","));

	std::copy(apY.begin(), apY.end() - 1,
		std::ostream_iterator<double>(oss4, ","));

	oss3 << _apY.back();
	oss4 << apY.back();

	std::string s3 = oss1.str(),
		s4 = oss2.str();


	REQUIRE(s3 == s4);

}
