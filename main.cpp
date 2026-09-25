#include <iostream>
#include "src/physics/DipoleField.hpp"
#include "src/math/Vector3.hpp"
#include "src/physics/Constants.hpp"

// TIP To <b>Run</b> code, press <shortcut actionId="Run"/> or click the <icon src="AllIcons.Actions.Execute"/> icon in the gutter.
int main() {
    Vector3 v1(C::earthRadius_m, 0, 0);
    Vector3 v2(0, 0, C::earthRadius_m);

    DipoleField DP;
    std::cout << DP.evaluate(v1) / C::b0_T << std::endl;
    std::cout << DP.evaluate(v2) / C::b0_T << std::endl;


    return 0;
    // TIP See CLion help at <a href="https://www.jetbrains.com/help/clion/">jetbrains.com/help/clion/</a>. Also, you can try interactive lessons for CLion by selecting 'Help | Learn IDE Features' from the main menu.
}