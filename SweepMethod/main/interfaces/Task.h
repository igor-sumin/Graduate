#pragma once

class Task {
protected:
    enum class TASK { NON_TASK, TASK_7 };

    static inline std::string getValue(TASK v) {
        switch (v) {
            case TASK::NON_TASK:
                return "not a task";
            case TASK::TASK_7:
                return "task 7";
            default:
                return "[Unknown TASK type]";
        }
    }

    virtual void defineDataByTask() {
        switch (task) {
            case TASK::TASK_7:
                this->defineDataByTask7();
                return;

            case TASK::NON_TASK:
                this->defineDataByNonTask();
                return;

            default:
                throw std::invalid_argument(Instrumental::getValue(task));
        }
    }

    virtual void defineDataByTask7() = 0;

    virtual void defineDataByNonTask() = 0;

    virtual void prepareData() = 0;

    virtual bool checkData() const = 0;
};