#pragma once
#include <SQLiteCpp/SQLiteCpp.h>
#include <chrono>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <boost/asio/ip/host_name.hpp>
#include <db.hpp>
#include <unistd.h>


namespace SQLite {

template <typename = void>
void create_database(const std::string& name) {
	SQLite::Database db(name, SQLite::OPEN_READWRITE|SQLite::OPEN_CREATE);
    SQLite::Transaction transaction(db);
    db.exec(DB::schema.c_str());
    transaction.commit();
}

template <typename T>
Data<T> get_from_database(const std::string& name, const std::string& shape, float radius, float height, const std::string& function, unsigned i, unsigned j) {
	try {
	SQLite::Database db(name, SQLite::OPEN_READWRITE);
	SQLite::Statement query(db, "SELECT value FROM Tuples WHERE shape=? AND radius=? AND height=? AND function=? AND i_value=? AND j_value=? AND value_type=?;");

	query.bind(1, shape);
	query.bind(2, radius);
	query.bind(3, height);
	query.bind(4, function);
	query.bind(5, i);
	query.bind(6, j);
	query.bind(7, typeid(T).name());

	bool got = query.executeStep();
	if (got == true) {
		std::string str = query.getColumn(0);
		return Data<T>(T(str), true);
	} else return Data<T>();
	} catch (std::exception& e) {
		// std::cout << "Exception found when reading. Returning empty Data<T>()" << std::endl;
		// std::cout << e.what() << std::endl;
		return Data<T>();
	}
}

template <typename T>
void save_at_database(const std::string& name, const std::string& shape, float radius, float height, const std::string& function, unsigned i, unsigned j, const T& result, int counter=0) {
	try {
		SQLite::Database db(name, SQLite::OPEN_READWRITE);
		SQLite::Transaction transaction(db);

		unsigned timestamp = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::system_clock::now().time_since_epoch()).count();
		std::string localhost = boost::asio::ip::host_name();

		std::stringstream sql;
		sql << "INSERT INTO Tuples (calculated_at, calculated_by, shape, radius, height, function, i_value, j_value, value_type, value) VALUES (";
		sql << timestamp << ", ";
		sql << "'" << localhost << "', ";
		sql << "'" << shape << "', ";
		sql << radius << ", ";
		sql << height << ", ";
		sql << "'" << function << "', ";
		sql << i << ", ";
		sql << j << ", ";
		sql << "'" << typeid(result).name() << "', ";
		sql << "'" << std::setprecision(std::numeric_limits<T>::max_digits10) << result << "'";
		sql << ");";

		db.exec(sql.str());
		transaction.commit();
	}

	// https://www.sqlite.org/rescode.html
	catch (SQLite::Exception& e) {
		/*std::cout << "Attempt: " << counter << std::endl;
		std::cout << "exception: " << e.what() << std::endl; 
		std::cout << "error code: " << e.getErrorCode() << std::endl; 
		std::cout << "extended code: " << e.getExtendedErrorCode() << std::endl; 
		std::cout << "str error: " << e.getErrorStr() << std::endl; */

		// SQLITE_BUSY
		int code = e.getErrorCode();
		if (code == 5) {
			usleep(100000);
			save_at_database(name, shape, radius, height, function, i, j, result, ++counter);

			if (counter > 500) throw e;
		}

		else {
			throw e;
		}
	}
}


} // SQLite Namespace.
