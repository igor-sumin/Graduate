#include <test/ParallelAlgorithmComponentTest.h>
#include <test/SerialAlgorithmComponentTest.h>
#include <test/common/MeasureTimeComponentTest.h>

int main(int argc, char** argv) {
	ParallelAlgorithmComponentTest p;
    p.execute();

//    SerialAlgorithmComponentTest s;
//    s.execute();
//
//    MeasureTimeComponentTest m;
//    m.execute();

	return 0;
}