#pragma once
#include <pqxx/pqxx>
#include <string>
#include <sstream>
#include <limits>
#include <iomanip>
#include <chrono>
#include <db.hpp>
#include <boost/asio/ip/host_name.hpp>

/*
Actually.. it is better to run all DB locally.
And, to have a third party program (python) syncronizing data with a PostgreSQL database.
*/


namespace PostGreSQL {
const std::string login_credentials = "dbname = fef user = postgres hostaddr = 127.0.0.1 port = 5432";

template <typename = void>
void create_database(const std::string& login) {
	pqxx::connection conn(login);
	pqxx::work transaction(conn);

	transaction.exec(DB::schema);
	transaction.commit();
	conn.disconnect();
}

template <typename T>
Data<T> get_from_database(const std::string& login, const std::string& shape, float radius, float height, const std::string& function, unsigned i, unsigned j) {
	try {
	pqxx::connection conn(login);
	pqxx::nontransaction obj(conn);

	std::stringstream sql;
	sql << "SELECT value FROM Tuples WHERE ";
	sql << "shape = " << obj.quote(shape) << " AND ";
	sql << "radius = " << radius << " AND ";
	sql << "height = " << height << " AND ";
	sql << "function = " << obj.quote(function) << " AND ";
	sql << "i_value = " << i << " AND ";
	sql << "j_value = " << j << " AND ";
	sql << "value_type = " << obj.quote(typeid(T).name()) << ";";

	pqxx::result data(obj.exec(sql.str()));
	if (data.empty()) {
		return Data<T>(T(), false);
	}
	
	std::string string_result = data[0][0].as<std::string>();
	conn.disconnect();
	
	return Data<T>(T(string_result), true);
	} catch(std::exception& e) {return Data<T>();}
}


template <typename T>
void save_at_database(const std::string& login, const std::string& shape, float radius, float height, const std::string& function, unsigned i, unsigned j, const T& result) {
	pqxx::connection conn(login);
	pqxx::work transaction(conn);

	unsigned timestamp = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
	std::string localhost = boost::asio::ip::host_name();

	std::stringstream sql;
	sql << "INSERT INTO Tuples (calculated_at, calculated_by, shape, radius, height, function, i_value, j_value, value_type, value) VALUES (";
	sql << timestamp << ", ";
	sql << transaction.quote(localhost) << ", ";
	sql << transaction.quote(shape) << ", ";
	sql << radius << ", ";
	sql << height << ", ";
	sql << transaction.quote(function) << ", ";
	sql << i << ", ";
	sql << j << ", ";
	sql << transaction.quote(typeid(result).name()) << ", ";
	sql << std::setprecision(std::numeric_limits<T>::max_digits10) << result;
	sql << ");";

	transaction.exec(sql.str());
	transaction.commit();
	conn.disconnect();
}

} // PostGreSQL namespace.
