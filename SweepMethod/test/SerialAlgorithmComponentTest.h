#pragma once

#include <functional>
#include <iostream>
#include <vector>
#include <test/common/BaseComponentTest.h>


class SerialAlgorithmComponentTest final : public BaseComponentTest {
private:
    double h;
    vec x;
    vec A, C, B;

    vec Phi;
    pairs kappa, mu, gamma;

    void prepareSerialDataForTest(const SerialSweepMethod& sweepMethod) {
        std::tie(v, u,
                 N, node, h, x,
                 A, C, B,
                 Phi, kappa, mu, gamma) = sweepMethod.getFields();
    }

    void setSerialFields(SerialSweepMethod& sweepMethod) {
        sweepMethod.setAllFields(v, u,
                                 N, node, h, x,
                                 A, C, B,
                                 Phi, kappa, mu, gamma);
    }

public:
    void testModelTask() {
        std::cout << "III) Module 7. Model task for estimating computational error:\n";

        // creation
        SerialSweepMethod ssm(5);
        this->prepareSerialDataForTest(ssm);

        for (size_t i = 0; i < node; i++) {
            u[i] = 10 + 90 * x[i] * x[i];
        }

        double total = 12. / (h * h);
        A.assign(N, total);
        C.assign(N, 2. * total + 5.);
        B.assign(N, total);

        mu = std::make_pair(10., 100.);
        kappa = std::make_pair(0., 0.);
        gamma = std::make_pair(1., 1.);

        for (size_t i = 0; i < node - 2; i++) {
            Phi[i] = 450. * x[i + 1] * x[i + 1] - 2110.;
        }
        this->setSerialFields(ssm);

        // execution
        v = ssm.run();
        ASSERT_FOR_DOUBLES(u, v);

        // check
        Phi.assign(node, 1);
        Phi[0] = mu.first; Phi[node - 1] = mu.second;
        for (size_t i = 1; i < node - 1; i++) {
            Phi[i] = 2110. - 450. * x[i] * x[i];
        }

        this->setSerialFields(ssm);
        vec res = Instrumental::calcMatrVecMult(ssm.createMatr(), v);

        printf("The scheme (SLAU) is solved with a discrepancy ||R|| = %f\n", ssm.calcR(res, Phi));
        printf("Estimation of the scheme error Z = %f\n", ssm.calcZ());
    }

    void execute() {
        std::vector<std::function<void()>> tests = {
            [this]() { this->testModelTask(); }
        };

        // BaseComponentTest::execute(tests);
    }
};