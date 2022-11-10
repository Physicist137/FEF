#pragma once
#include <string>

template <typename T>
struct Data {
	T result;
	bool found;

	explicit Data(const T& res = T(), bool f = false) : result(res), found(f) {}
};


namespace DB {


// https://stackoverflow.com/questions/787722/postgresql-autoincrement
const std::string schema = "CREATE TABLE Tuples(" \
	"id integer not null primary key autoincrement,"
	"calculated_at integer not null, " \
	"calculated_by text not null, " \
	"shape text not null, " \
	"radius real not null, " \
	"height real not null, " \
	"function text not null, " \
	"i_value integer, " \
	"j_value integer, " \
	"value_type text not null, " \
	"value text not null);";

} // DB Namespace.
