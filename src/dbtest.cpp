#include <postgresqldb.hpp>
#include <sqlitedb.hpp>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/multiprecision/float128.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <chrono>

void Sqlite() {
	using Float = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<8192>>;
	using SmallFloat = boost::multiprecision::float128;
	const std::string login = "tester.db";
	const std::string shape = "prolate_spheroid";
	SQLite::save_at_database<SmallFloat>(login, shape, 1.0, 2.0, "test", 2, 3, 3.1415926535);

	auto a = SQLite::get_from_database<SmallFloat>(login, shape, 1.0, 2.0, "test", 2, 3);
	std::cout << a.found << ") " << a.result << std::endl;
}

void Postgres() {
	using Float = boost::multiprecision::number<boost::multiprecision::cpp_bin_float<8192>>;
	using SmallFloat = boost::multiprecision::float128;

	const std::string login = "dbname = fef user = postgres hostaddr = 127.0.0.1 port = 5432";
	const std::string shape = "prolate_spheroid";

	PostGreSQL::create_database(login);
	PostGreSQL::save_at_database<SmallFloat>(login, shape, 1.0, 2.0, "test", 4, 0, SmallFloat("45.498167165784"));

	auto a = PostGreSQL::get_from_database<SmallFloat>(login, shape, 1.0, 2.0, "test", 4, 0);

	std::cout << a.found << ") " << a.result << std::endl;
}

int main() {
	Sqlite();
}
