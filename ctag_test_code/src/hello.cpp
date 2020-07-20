//#include "../interface/hello.C"
#include "../interface/ctag_reshaping.C"

int main()
{
    printf("Hello World!\n");

    retrieve_scale_factor sf;
    TRandom3 random(1234);

    double cvsl = random.Rndm();
    double cvsb = random.Rndm();
    TString type_flavour = "b";
    TString name = "SF" + type_flavour + "_hist";

    //sf.debug_mode();
    double scale_factor = sf.get_scale_factor(name, cvsl, cvsb);

    std::cout << "This is the end of file" << std::endl;
    return 0;
}
