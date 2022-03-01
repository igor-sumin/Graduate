#include <test/ParallelAlgorithmComponentTest.h>
#include <test/SerialAlgorithmComponentTest.h>
#include <test/common/InstrumentalComponentTest.h>

int main(int argc, char** argv) {
	ParallelAlgorithmComponentTest p;
    p.execute();

    SerialAlgorithmComponentTest s;
    s.execute();

    InstrumentalComponentTest i;
    i.execute();

	return 0;
}