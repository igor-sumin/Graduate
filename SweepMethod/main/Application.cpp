
#include <test/ParallelAlgorithmComponentTest.h>
#include <test/SerialAlgorithmComponentTest.h>

int main(int argc, char** argv) {
	ParallelAlgorithmComponentTest p;
    p.execute();

    SerialAlgorithmComponentTest s;
    // s.execute();

	return 0;
}