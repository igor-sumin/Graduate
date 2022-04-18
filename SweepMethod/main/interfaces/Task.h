#pragma once

class Task {
public:
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

protected:
    virtual void defineDataByTask7() = 0;

    virtual void defineDataByNonTask() = 0;

    virtual void prepareData() = 0;

    virtual bool checkData() const = 0;
};