#include <utility>
#include <fstream>
#include <iomanip>

#pragma once

class FWork {
private:
    str fullPath;

protected:
    bool write(std::ofstream& output, const vec2d<double>& u, const str& step) {
        if (!output.is_open()) {
            return false;
        }

        output << std::setprecision(5);
        output << "Шаг s -------- " << step << "\n";

        loop3(k) {
            output << "компонента --- " << k << "\n";
            // работа с 3 компонентами
            for (size_t j = 0; j < u[k][0].size(); j++) {
                for (size_t i = 0; i < u[k].size(); i++) {
                    output << std::setw(9) << u[k][i][j] << " ";
                } output << "\n";
            } output << "\n";
        } output << "\n";

        output.close();
        return true;
    }

public:
    explicit FWork() {
        fullPath = AppConstansts::FULL_PATH;
    }

    // Запись данных в файл
    bool fwrite(const vec2d<double>& uPhase, size_t s, const str& layer) {
        std::ofstream output(fullPath + layer + ".txt", std::ios::app);

        str step = (layer == AppConstansts::MAIN_LAYER)
                    ? std::to_string(s)
                    : std::to_string(s) + ".5";

        this->write(output, uPhase, step);

        return true;
    }

    // Чтение данных из файла paths
    bool fread(const str& layer) {
        std::ifstream input(fullPath + layer + ".txt");

        if (!input.is_open()) {
            return false;
        }

        str line;
        while (std::getline(input, line)) {
            std::cout << line << "\n";
        }

        input.close();
        return true;
    }
};