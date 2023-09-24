#include "table.h"

void printTable(const std::vector<std::string>& headers, const std::vector<std::vector<std::string>>& data) {
    if (headers.empty() || data.empty()) {
        std::cout << "Table is empty." << std::endl;
        return;
    }

    size_t numColumns = headers.size();
    std::vector<size_t> columnWidths(numColumns);
    for (size_t i = 0; i < numColumns; ++i) {
        columnWidths[i] = headers[i].length();
        for (const auto& row : data) {
            if (row.size() > i && row[i].length() > columnWidths[i]) {
                columnWidths[i] = row[i].length();
            }
        }
    }

    for (size_t i = 0; i < numColumns; ++i) {
        std::cout << std::setw(columnWidths[i]) << headers[i] << " | ";
    }
    std::cout << std::endl;

    for (size_t i = 0; i < numColumns; ++i) {
        std::cout << std::string(columnWidths[i], '-') << "-+-";
    }
    std::cout << std::endl;

    for (const auto& row : data) {
        if (row.size() != numColumns) {
            std::cout << "Invalid row size." << std::endl;
            return;
        }
        for (size_t i = 0; i < numColumns; ++i) {
            std::cout << std::setw(columnWidths[i]) << row[i] << " | ";
        }
        std::cout << std::endl;
    }
}

