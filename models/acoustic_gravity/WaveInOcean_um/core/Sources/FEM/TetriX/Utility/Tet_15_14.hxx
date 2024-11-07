#pragma once


// -----------------------------------------------------------------------------------------//
namespace OndoMathX
{

    namespace TetriX
    {
 
        namespace Tet_15_14
        {
        
            inline void TransposeGradientInterpolation(const Real * in, Real *out)
            {
                
                const Real c1 = TetriX::Tet_14::_c1;
                const Real c2 = TetriX::Tet_14::_c2;
                const Real d = TetriX::Tet_14::_d;

                const Real c1_2 = c1*c1;
                const Real c2_2 = c2*c2;
                const Real d_2 = d*d;

                const Real c1_3 = c1*c1*c1;
                const Real c2_3 = c2*c2*c2;
                const Real d_3 = d*d*d;

                Real Buffer;

                for (Index i=0;i<15;++i) 
                {
                    out[i]=0.0;
                }
Buffer = (-1.5846000073)*in[0]; 
out[0] += Buffer; 
 
Buffer = (-0.2748830499)*in[0]; 
out[1] += Buffer; 
 
Buffer = (0.1533686378)*in[0]; 
out[2] += Buffer; 
 
Buffer = (0.1533686378)*in[0]; 
out[3] += Buffer; 
 
Buffer = (1.2892868929)*in[0]; 
out[4] += Buffer; 
 
Buffer = (-0.2591746961)*in[0]; 
out[5] += Buffer; 
 
Buffer = (-0.7946608590)*in[0]; 
out[6] += Buffer; 
 
Buffer = (-0.7946608590)*in[0]; 
out[7] += Buffer; 
 
Buffer = (-0.2591746961)*in[0]; 
out[8] += Buffer; 
 
Buffer = (0.1731135469)*in[0]; 
out[9] += Buffer; 
 
Buffer = (0.9908122592)*in[0]; 
out[10] += Buffer; 
 
Buffer = (0.9908122592)*in[0]; 
out[11] += Buffer; 
 
Buffer = (-0.8164535405)*in[0]; 
out[12] += Buffer; 
 
Buffer = (-0.3520629014)*in[0]; 
out[13] += Buffer; 
 
Buffer = (1.3849083756)*in[0]; 
out[14] += Buffer; 
 
Buffer = (-1.5846000073)*in[1]; 
out[0] += Buffer; 
 
Buffer = (0.1533686378)*in[1]; 
out[1] += Buffer; 
 
Buffer = (-0.2748830499)*in[1]; 
out[2] += Buffer; 
 
Buffer = (0.1533686378)*in[1]; 
out[3] += Buffer; 
 
Buffer = (-0.7946608590)*in[1]; 
out[4] += Buffer; 
 
Buffer = (-0.2591746961)*in[1]; 
out[5] += Buffer; 
 
Buffer = (1.2892868929)*in[1]; 
out[6] += Buffer; 
 
Buffer = (-0.7946608590)*in[1]; 
out[7] += Buffer; 
 
Buffer = (0.1731135469)*in[1]; 
out[8] += Buffer; 
 
Buffer = (-0.2591746961)*in[1]; 
out[9] += Buffer; 
 
Buffer = (0.9908122592)*in[1]; 
out[10] += Buffer; 
 
Buffer = (-0.8164535405)*in[1]; 
out[11] += Buffer; 
 
Buffer = (0.9908122592)*in[1]; 
out[12] += Buffer; 
 
Buffer = (-0.3520629014)*in[1]; 
out[13] += Buffer; 
 
Buffer = (1.3849083756)*in[1]; 
out[14] += Buffer; 
 
Buffer = (-1.5846000073)*in[2]; 
out[0] += Buffer; 
 
Buffer = (0.1533686378)*in[2]; 
out[1] += Buffer; 
 
Buffer = (0.1533686378)*in[2]; 
out[2] += Buffer; 
 
Buffer = (-0.2748830499)*in[2]; 
out[3] += Buffer; 
 
Buffer = (-0.7946608590)*in[2]; 
out[4] += Buffer; 
 
Buffer = (0.1731135469)*in[2]; 
out[5] += Buffer; 
 
Buffer = (-0.7946608590)*in[2]; 
out[6] += Buffer; 
 
Buffer = (1.2892868929)*in[2]; 
out[7] += Buffer; 
 
Buffer = (-0.2591746961)*in[2]; 
out[8] += Buffer; 
 
Buffer = (-0.2591746961)*in[2]; 
out[9] += Buffer; 
 
Buffer = (-0.8164535405)*in[2]; 
out[10] += Buffer; 
 
Buffer = (0.9908122592)*in[2]; 
out[11] += Buffer; 
 
Buffer = (0.9908122592)*in[2]; 
out[12] += Buffer; 
 
Buffer = (-0.3520629014)*in[2]; 
out[13] += Buffer; 
 
Buffer = (1.3849083756)*in[2]; 
out[14] += Buffer; 
 
Buffer = (0.2748830499)*in[3]; 
out[0] += Buffer; 
 
Buffer = (1.5846000073)*in[3]; 
out[1] += Buffer; 
 
Buffer = (-0.1533686378)*in[3]; 
out[2] += Buffer; 
 
Buffer = (-0.1533686378)*in[3]; 
out[3] += Buffer; 
 
Buffer = (-1.2892868929)*in[3]; 
out[4] += Buffer; 
 
Buffer = (0.7946608590)*in[3]; 
out[5] += Buffer; 
 
Buffer = (0.2591746961)*in[3]; 
out[6] += Buffer; 
 
Buffer = (0.2591746961)*in[3]; 
out[7] += Buffer; 
 
Buffer = (0.7946608590)*in[3]; 
out[8] += Buffer; 
 
Buffer = (-0.1731135469)*in[3]; 
out[9] += Buffer; 
 
Buffer = (-0.9908122592)*in[3]; 
out[10] += Buffer; 
 
Buffer = (-0.9908122592)*in[3]; 
out[11] += Buffer; 
 
Buffer = (0.3520629014)*in[3]; 
out[12] += Buffer; 
 
Buffer = (0.8164535405)*in[3]; 
out[13] += Buffer; 
 
Buffer = (-1.3849083756)*in[3]; 
out[14] += Buffer; 
 
Buffer = (0.4282516877)*in[4]; 
out[0] += Buffer; 
 
Buffer = (-0.4282516877)*in[4]; 
out[2] += Buffer; 
 
Buffer = (-2.0839477519)*in[4]; 
out[4] += Buffer; 
 
Buffer = (2.0839477519)*in[4]; 
out[5] += Buffer; 
 
Buffer = (0.4322882431)*in[4]; 
out[7] += Buffer; 
 
Buffer = (-0.4322882431)*in[4]; 
out[9] += Buffer; 
 
Buffer = (-1.8072657997)*in[4]; 
out[11] += Buffer; 
 
Buffer = (1.8072657997)*in[4]; 
out[13] += Buffer; 
 
Buffer = (0.4282516877)*in[5]; 
out[0] += Buffer; 
 
Buffer = (-0.4282516877)*in[5]; 
out[3] += Buffer; 
 
Buffer = (-2.0839477519)*in[5]; 
out[4] += Buffer; 
 
Buffer = (0.4322882431)*in[5]; 
out[6] += Buffer; 
 
Buffer = (2.0839477519)*in[5]; 
out[8] += Buffer; 
 
Buffer = (-0.4322882431)*in[5]; 
out[9] += Buffer; 
 
Buffer = (-1.8072657997)*in[5]; 
out[10] += Buffer; 
 
Buffer = (1.8072657997)*in[5]; 
out[13] += Buffer; 
 
Buffer = (0.4282516877)*in[6]; 
out[0] += Buffer; 
 
Buffer = (-0.4282516877)*in[6]; 
out[1] += Buffer; 
 
Buffer = (2.0839477519)*in[6]; 
out[5] += Buffer; 
 
Buffer = (-2.0839477519)*in[6]; 
out[6] += Buffer; 
 
Buffer = (0.4322882431)*in[6]; 
out[7] += Buffer; 
 
Buffer = (-0.4322882431)*in[6]; 
out[8] += Buffer; 
 
Buffer = (-1.8072657997)*in[6]; 
out[12] += Buffer; 
 
Buffer = (1.8072657997)*in[6]; 
out[13] += Buffer; 
 
Buffer = (0.2748830499)*in[7]; 
out[0] += Buffer; 
 
Buffer = (-0.1533686378)*in[7]; 
out[1] += Buffer; 
 
Buffer = (1.5846000073)*in[7]; 
out[2] += Buffer; 
 
Buffer = (-0.1533686378)*in[7]; 
out[3] += Buffer; 
 
Buffer = (0.2591746961)*in[7]; 
out[4] += Buffer; 
 
Buffer = (0.7946608590)*in[7]; 
out[5] += Buffer; 
 
Buffer = (-1.2892868929)*in[7]; 
out[6] += Buffer; 
 
Buffer = (0.2591746961)*in[7]; 
out[7] += Buffer; 
 
Buffer = (-0.1731135469)*in[7]; 
out[8] += Buffer; 
 
Buffer = (0.7946608590)*in[7]; 
out[9] += Buffer; 
 
Buffer = (-0.9908122592)*in[7]; 
out[10] += Buffer; 
 
Buffer = (0.3520629014)*in[7]; 
out[11] += Buffer; 
 
Buffer = (-0.9908122592)*in[7]; 
out[12] += Buffer; 
 
Buffer = (0.8164535405)*in[7]; 
out[13] += Buffer; 
 
Buffer = (-1.3849083756)*in[7]; 
out[14] += Buffer; 
 
Buffer = (0.4282516877)*in[8]; 
out[0] += Buffer; 
 
Buffer = (-0.4282516877)*in[8]; 
out[3] += Buffer; 
 
Buffer = (0.4322882431)*in[8]; 
out[4] += Buffer; 
 
Buffer = (-2.0839477519)*in[8]; 
out[6] += Buffer; 
 
Buffer = (-0.4322882431)*in[8]; 
out[8] += Buffer; 
 
Buffer = (2.0839477519)*in[8]; 
out[9] += Buffer; 
 
Buffer = (-1.8072657997)*in[8]; 
out[10] += Buffer; 
 
Buffer = (1.8072657997)*in[8]; 
out[13] += Buffer; 
 
Buffer = (0.4282516877)*in[9]; 
out[0] += Buffer; 
 
Buffer = (-0.4282516877)*in[9]; 
out[1] += Buffer; 
 
Buffer = (-0.4322882431)*in[9]; 
out[5] += Buffer; 
 
Buffer = (0.4322882431)*in[9]; 
out[6] += Buffer; 
 
Buffer = (-2.0839477519)*in[9]; 
out[7] += Buffer; 
 
Buffer = (2.0839477519)*in[9]; 
out[8] += Buffer; 
 
Buffer = (-1.8072657997)*in[9]; 
out[12] += Buffer; 
 
Buffer = (1.8072657997)*in[9]; 
out[13] += Buffer; 
 
Buffer = (0.4282516877)*in[10]; 
out[0] += Buffer; 
 
Buffer = (-0.4282516877)*in[10]; 
out[2] += Buffer; 
 
Buffer = (0.4322882431)*in[10]; 
out[4] += Buffer; 
 
Buffer = (-0.4322882431)*in[10]; 
out[5] += Buffer; 
 
Buffer = (-2.0839477519)*in[10]; 
out[7] += Buffer; 
 
Buffer = (2.0839477519)*in[10]; 
out[9] += Buffer; 
 
Buffer = (-1.8072657997)*in[10]; 
out[11] += Buffer; 
 
Buffer = (1.8072657997)*in[10]; 
out[13] += Buffer; 
 
Buffer = (0.2748830499)*in[11]; 
out[0] += Buffer; 
 
Buffer = (-0.1533686378)*in[11]; 
out[1] += Buffer; 
 
Buffer = (-0.1533686378)*in[11]; 
out[2] += Buffer; 
 
Buffer = (1.5846000073)*in[11]; 
out[3] += Buffer; 
 
Buffer = (0.2591746961)*in[11]; 
out[4] += Buffer; 
 
Buffer = (-0.1731135469)*in[11]; 
out[5] += Buffer; 
 
Buffer = (0.2591746961)*in[11]; 
out[6] += Buffer; 
 
Buffer = (-1.2892868929)*in[11]; 
out[7] += Buffer; 
 
Buffer = (0.7946608590)*in[11]; 
out[8] += Buffer; 
 
Buffer = (0.7946608590)*in[11]; 
out[9] += Buffer; 
 
Buffer = (0.3520629014)*in[11]; 
out[10] += Buffer; 
 
Buffer = (-0.9908122592)*in[11]; 
out[11] += Buffer; 
 
Buffer = (-0.9908122592)*in[11]; 
out[12] += Buffer; 
 
Buffer = (0.8164535405)*in[11]; 
out[13] += Buffer; 
 
Buffer = (-1.3849083756)*in[11]; 
out[14] += Buffer; 
 
Buffer = (0.0805491062)*in[12]; 
out[0] += Buffer; 
 
Buffer = (0.1733620808)*in[12]; 
out[1] += Buffer; 
 
Buffer = (-0.1329888608)*in[12]; 
out[2] += Buffer; 
 
Buffer = (-0.1329888608)*in[12]; 
out[3] += Buffer; 
 
Buffer = (0.0897361785)*in[12]; 
out[4] += Buffer; 
 
Buffer = (0.2390823071)*in[12]; 
out[5] += Buffer; 
 
Buffer = (0.0715962681)*in[12]; 
out[6] += Buffer; 
 
Buffer = (0.0715962681)*in[12]; 
out[7] += Buffer; 
 
Buffer = (0.2390823071)*in[12]; 
out[8] += Buffer; 
 
Buffer = (-0.7532323114)*in[12]; 
out[9] += Buffer; 
 
Buffer = (0.4978729532)*in[12]; 
out[10] += Buffer; 
 
Buffer = (0.4978729532)*in[12]; 
out[11] += Buffer; 
 
Buffer = (-0.0673924284)*in[12]; 
out[12] += Buffer; 
 
Buffer = (5.1517105306)*in[12]; 
out[13] += Buffer; 
 
Buffer = (-6.0258584915)*in[12]; 
out[14] += Buffer; 
 
Buffer = (0.0805491062)*in[13]; 
out[0] += Buffer; 
 
Buffer = (-0.1329888608)*in[13]; 
out[1] += Buffer; 
 
Buffer = (0.1733620808)*in[13]; 
out[2] += Buffer; 
 
Buffer = (-0.1329888608)*in[13]; 
out[3] += Buffer; 
 
Buffer = (0.0715962681)*in[13]; 
out[4] += Buffer; 
 
Buffer = (0.2390823071)*in[13]; 
out[5] += Buffer; 
 
Buffer = (0.0897361785)*in[13]; 
out[6] += Buffer; 
 
Buffer = (0.0715962681)*in[13]; 
out[7] += Buffer; 
 
Buffer = (-0.7532323114)*in[13]; 
out[8] += Buffer; 
 
Buffer = (0.2390823071)*in[13]; 
out[9] += Buffer; 
 
Buffer = (0.4978729532)*in[13]; 
out[10] += Buffer; 
 
Buffer = (-0.0673924284)*in[13]; 
out[11] += Buffer; 
 
Buffer = (0.4978729532)*in[13]; 
out[12] += Buffer; 
 
Buffer = (5.1517105306)*in[13]; 
out[13] += Buffer; 
 
Buffer = (-6.0258584915)*in[13]; 
out[14] += Buffer; 
 
Buffer = (0.0805491062)*in[14]; 
out[0] += Buffer; 
 
Buffer = (-0.1329888608)*in[14]; 
out[1] += Buffer; 
 
Buffer = (-0.1329888608)*in[14]; 
out[2] += Buffer; 
 
Buffer = (0.1733620808)*in[14]; 
out[3] += Buffer; 
 
Buffer = (0.0715962681)*in[14]; 
out[4] += Buffer; 
 
Buffer = (-0.7532323114)*in[14]; 
out[5] += Buffer; 
 
Buffer = (0.0715962681)*in[14]; 
out[6] += Buffer; 
 
Buffer = (0.0897361785)*in[14]; 
out[7] += Buffer; 
 
Buffer = (0.2390823071)*in[14]; 
out[8] += Buffer; 
 
Buffer = (0.2390823071)*in[14]; 
out[9] += Buffer; 
 
Buffer = (-0.0673924284)*in[14]; 
out[10] += Buffer; 
 
Buffer = (0.4978729532)*in[14]; 
out[11] += Buffer; 
 
Buffer = (0.4978729532)*in[14]; 
out[12] += Buffer; 
 
Buffer = (5.1517105306)*in[14]; 
out[13] += Buffer; 
 
Buffer = (-6.0258584915)*in[14]; 
out[14] += Buffer; 
 
Buffer = (-0.1733620808)*in[15]; 
out[0] += Buffer; 
 
Buffer = (-0.0805491062)*in[15]; 
out[1] += Buffer; 
 
Buffer = (0.1329888608)*in[15]; 
out[2] += Buffer; 
 
Buffer = (0.1329888608)*in[15]; 
out[3] += Buffer; 
 
Buffer = (-0.0897361785)*in[15]; 
out[4] += Buffer; 
 
Buffer = (-0.0715962681)*in[15]; 
out[5] += Buffer; 
 
Buffer = (-0.2390823071)*in[15]; 
out[6] += Buffer; 
 
Buffer = (-0.2390823071)*in[15]; 
out[7] += Buffer; 
 
Buffer = (-0.0715962681)*in[15]; 
out[8] += Buffer; 
 
Buffer = (0.7532323114)*in[15]; 
out[9] += Buffer; 
 
Buffer = (-0.4978729532)*in[15]; 
out[10] += Buffer; 
 
Buffer = (-0.4978729532)*in[15]; 
out[11] += Buffer; 
 
Buffer = (-5.1517105306)*in[15]; 
out[12] += Buffer; 
 
Buffer = (0.0673924284)*in[15]; 
out[13] += Buffer; 
 
Buffer = (6.0258584915)*in[15]; 
out[14] += Buffer; 
 
Buffer = (-0.3063509417)*in[16]; 
out[0] += Buffer; 
 
Buffer = (0.3063509417)*in[16]; 
out[2] += Buffer; 
 
Buffer = (-0.0181399103)*in[16]; 
out[4] += Buffer; 
 
Buffer = (0.0181399103)*in[16]; 
out[5] += Buffer; 
 
Buffer = (-0.9923146186)*in[16]; 
out[7] += Buffer; 
 
Buffer = (0.9923146186)*in[16]; 
out[9] += Buffer; 
 
Buffer = (-0.5652653816)*in[16]; 
out[11] += Buffer; 
 
Buffer = (0.5652653816)*in[16]; 
out[13] += Buffer; 
 
Buffer = (-0.3063509417)*in[17]; 
out[0] += Buffer; 
 
Buffer = (0.3063509417)*in[17]; 
out[3] += Buffer; 
 
Buffer = (-0.0181399103)*in[17]; 
out[4] += Buffer; 
 
Buffer = (-0.9923146186)*in[17]; 
out[6] += Buffer; 
 
Buffer = (0.0181399103)*in[17]; 
out[8] += Buffer; 
 
Buffer = (0.9923146186)*in[17]; 
out[9] += Buffer; 
 
Buffer = (-0.5652653816)*in[17]; 
out[10] += Buffer; 
 
Buffer = (0.5652653816)*in[17]; 
out[13] += Buffer; 
 
Buffer = (-0.3063509417)*in[18]; 
out[0] += Buffer; 
 
Buffer = (0.3063509417)*in[18]; 
out[1] += Buffer; 
 
Buffer = (0.0181399103)*in[18]; 
out[5] += Buffer; 
 
Buffer = (-0.0181399103)*in[18]; 
out[6] += Buffer; 
 
Buffer = (-0.9923146186)*in[18]; 
out[7] += Buffer; 
 
Buffer = (0.9923146186)*in[18]; 
out[8] += Buffer; 
 
Buffer = (-0.5652653816)*in[18]; 
out[12] += Buffer; 
 
Buffer = (0.5652653816)*in[18]; 
out[13] += Buffer; 
 
Buffer = (-0.1733620808)*in[19]; 
out[0] += Buffer; 
 
Buffer = (0.1329888608)*in[19]; 
out[1] += Buffer; 
 
Buffer = (-0.0805491062)*in[19]; 
out[2] += Buffer; 
 
Buffer = (0.1329888608)*in[19]; 
out[3] += Buffer; 
 
Buffer = (-0.2390823071)*in[19]; 
out[4] += Buffer; 
 
Buffer = (-0.0715962681)*in[19]; 
out[5] += Buffer; 
 
Buffer = (-0.0897361785)*in[19]; 
out[6] += Buffer; 
 
Buffer = (-0.2390823071)*in[19]; 
out[7] += Buffer; 
 
Buffer = (0.7532323114)*in[19]; 
out[8] += Buffer; 
 
Buffer = (-0.0715962681)*in[19]; 
out[9] += Buffer; 
 
Buffer = (-0.4978729532)*in[19]; 
out[10] += Buffer; 
 
Buffer = (-5.1517105306)*in[19]; 
out[11] += Buffer; 
 
Buffer = (-0.4978729532)*in[19]; 
out[12] += Buffer; 
 
Buffer = (0.0673924284)*in[19]; 
out[13] += Buffer; 
 
Buffer = (6.0258584915)*in[19]; 
out[14] += Buffer; 
 
Buffer = (-0.3063509417)*in[20]; 
out[0] += Buffer; 
 
Buffer = (0.3063509417)*in[20]; 
out[3] += Buffer; 
 
Buffer = (-0.9923146186)*in[20]; 
out[4] += Buffer; 
 
Buffer = (-0.0181399103)*in[20]; 
out[6] += Buffer; 
 
Buffer = (0.9923146186)*in[20]; 
out[8] += Buffer; 
 
Buffer = (0.0181399103)*in[20]; 
out[9] += Buffer; 
 
Buffer = (-0.5652653816)*in[20]; 
out[10] += Buffer; 
 
Buffer = (0.5652653816)*in[20]; 
out[13] += Buffer; 
 
Buffer = (-0.3063509417)*in[21]; 
out[0] += Buffer; 
 
Buffer = (0.3063509417)*in[21]; 
out[1] += Buffer; 
 
Buffer = (0.9923146186)*in[21]; 
out[5] += Buffer; 
 
Buffer = (-0.9923146186)*in[21]; 
out[6] += Buffer; 
 
Buffer = (-0.0181399103)*in[21]; 
out[7] += Buffer; 
 
Buffer = (0.0181399103)*in[21]; 
out[8] += Buffer; 
 
Buffer = (-0.5652653816)*in[21]; 
out[12] += Buffer; 
 
Buffer = (0.5652653816)*in[21]; 
out[13] += Buffer; 
 
Buffer = (-0.3063509417)*in[22]; 
out[0] += Buffer; 
 
Buffer = (0.3063509417)*in[22]; 
out[2] += Buffer; 
 
Buffer = (-0.9923146186)*in[22]; 
out[4] += Buffer; 
 
Buffer = (0.9923146186)*in[22]; 
out[5] += Buffer; 
 
Buffer = (-0.0181399103)*in[22]; 
out[7] += Buffer; 
 
Buffer = (0.0181399103)*in[22]; 
out[9] += Buffer; 
 
Buffer = (-0.5652653816)*in[22]; 
out[11] += Buffer; 
 
Buffer = (0.5652653816)*in[22]; 
out[13] += Buffer; 
 
Buffer = (-0.1733620808)*in[23]; 
out[0] += Buffer; 
 
Buffer = (0.1329888608)*in[23]; 
out[1] += Buffer; 
 
Buffer = (0.1329888608)*in[23]; 
out[2] += Buffer; 
 
Buffer = (-0.0805491062)*in[23]; 
out[3] += Buffer; 
 
Buffer = (-0.2390823071)*in[23]; 
out[4] += Buffer; 
 
Buffer = (0.7532323114)*in[23]; 
out[5] += Buffer; 
 
Buffer = (-0.2390823071)*in[23]; 
out[6] += Buffer; 
 
Buffer = (-0.0897361785)*in[23]; 
out[7] += Buffer; 
 
Buffer = (-0.0715962681)*in[23]; 
out[8] += Buffer; 
 
Buffer = (-0.0715962681)*in[23]; 
out[9] += Buffer; 
 
Buffer = (-5.1517105306)*in[23]; 
out[10] += Buffer; 
 
Buffer = (-0.4978729532)*in[23]; 
out[11] += Buffer; 
 
Buffer = (-0.4978729532)*in[23]; 
out[12] += Buffer; 
 
Buffer = (0.0673924284)*in[23]; 
out[13] += Buffer; 
 
Buffer = (6.0258584915)*in[23]; 
out[14] += Buffer; 
 
Buffer = (-0.8241969448)*in[24]; 
out[0] += Buffer; 
 
Buffer = (0.8241969448)*in[24]; 
out[1] += Buffer; 
 
Buffer = (0.1571677714)*in[24]; 
out[5] += Buffer; 
 
Buffer = (-0.1571677714)*in[24]; 
out[6] += Buffer; 
 
Buffer = (-0.1571677714)*in[24]; 
out[7] += Buffer; 
 
Buffer = (0.1571677714)*in[24]; 
out[8] += Buffer; 
 
Buffer = (-0.0559058514)*in[24]; 
out[12] += Buffer; 
 
Buffer = (0.0559058514)*in[24]; 
out[13] += Buffer; 
 
Buffer = (-0.3003740274)*in[25]; 
out[0] += Buffer; 
 
Buffer = (0.5238229173)*in[25]; 
out[1] += Buffer; 
 
Buffer = (-0.1762864376)*in[25]; 
out[2] += Buffer; 
 
Buffer = (0.0219980970)*in[25]; 
out[3] += Buffer; 
 
Buffer = (-3.5297659265)*in[25]; 
out[4] += Buffer; 
 
Buffer = (-0.3901459189)*in[25]; 
out[5] += Buffer; 
 
Buffer = (-0.5473136903)*in[25]; 
out[6] += Buffer; 
 
Buffer = (0.1135037218)*in[25]; 
out[7] += Buffer; 
 
Buffer = (0.2706714932)*in[25]; 
out[8] += Buffer; 
 
Buffer = (-0.0188170046)*in[25]; 
out[9] += Buffer; 
 
Buffer = (4.1053953962)*in[25]; 
out[10] += Buffer; 
 
Buffer = (-1.4719104438)*in[25]; 
out[11] += Buffer; 
 
Buffer = (-0.4110279866)*in[25]; 
out[12] += Buffer; 
 
Buffer = (-0.3551221352)*in[25]; 
out[13] += Buffer; 
 
Buffer = (2.1653719455)*in[25]; 
out[14] += Buffer; 
 
Buffer = (-0.3003740274)*in[26]; 
out[0] += Buffer; 
 
Buffer = (0.5238229173)*in[26]; 
out[1] += Buffer; 
 
Buffer = (0.0219980970)*in[26]; 
out[2] += Buffer; 
 
Buffer = (-0.1762864376)*in[26]; 
out[3] += Buffer; 
 
Buffer = (-3.5297659265)*in[26]; 
out[4] += Buffer; 
 
Buffer = (0.2706714932)*in[26]; 
out[5] += Buffer; 
 
Buffer = (0.1135037218)*in[26]; 
out[6] += Buffer; 
 
Buffer = (-0.5473136903)*in[26]; 
out[7] += Buffer; 
 
Buffer = (-0.3901459189)*in[26]; 
out[8] += Buffer; 
 
Buffer = (-0.0188170046)*in[26]; 
out[9] += Buffer; 
 
Buffer = (-1.4719104438)*in[26]; 
out[10] += Buffer; 
 
Buffer = (4.1053953962)*in[26]; 
out[11] += Buffer; 
 
Buffer = (-0.4110279866)*in[26]; 
out[12] += Buffer; 
 
Buffer = (-0.3551221352)*in[26]; 
out[13] += Buffer; 
 
Buffer = (2.1653719455)*in[26]; 
out[14] += Buffer; 
 
Buffer = (-0.3003740274)*in[27]; 
out[0] += Buffer; 
 
Buffer = (-0.1762864376)*in[27]; 
out[1] += Buffer; 
 
Buffer = (0.5238229173)*in[27]; 
out[2] += Buffer; 
 
Buffer = (0.0219980970)*in[27]; 
out[3] += Buffer; 
 
Buffer = (-0.5473136903)*in[27]; 
out[4] += Buffer; 
 
Buffer = (-0.3901459189)*in[27]; 
out[5] += Buffer; 
 
Buffer = (-3.5297659265)*in[27]; 
out[6] += Buffer; 
 
Buffer = (0.1135037218)*in[27]; 
out[7] += Buffer; 
 
Buffer = (-0.0188170046)*in[27]; 
out[8] += Buffer; 
 
Buffer = (0.2706714932)*in[27]; 
out[9] += Buffer; 
 
Buffer = (4.1053953962)*in[27]; 
out[10] += Buffer; 
 
Buffer = (-0.4110279866)*in[27]; 
out[11] += Buffer; 
 
Buffer = (-1.4719104438)*in[27]; 
out[12] += Buffer; 
 
Buffer = (-0.3551221352)*in[27]; 
out[13] += Buffer; 
 
Buffer = (2.1653719455)*in[27]; 
out[14] += Buffer; 
 
Buffer = (-0.8241969448)*in[28]; 
out[0] += Buffer; 
 
Buffer = (0.8241969448)*in[28]; 
out[2] += Buffer; 
 
Buffer = (-0.1571677714)*in[28]; 
out[4] += Buffer; 
 
Buffer = (0.1571677714)*in[28]; 
out[5] += Buffer; 
 
Buffer = (-0.1571677714)*in[28]; 
out[7] += Buffer; 
 
Buffer = (0.1571677714)*in[28]; 
out[9] += Buffer; 
 
Buffer = (-0.0559058514)*in[28]; 
out[11] += Buffer; 
 
Buffer = (0.0559058514)*in[28]; 
out[13] += Buffer; 
 
Buffer = (-0.3003740274)*in[29]; 
out[0] += Buffer; 
 
Buffer = (0.0219980970)*in[29]; 
out[1] += Buffer; 
 
Buffer = (0.5238229173)*in[29]; 
out[2] += Buffer; 
 
Buffer = (-0.1762864376)*in[29]; 
out[3] += Buffer; 
 
Buffer = (0.1135037218)*in[29]; 
out[4] += Buffer; 
 
Buffer = (0.2706714932)*in[29]; 
out[5] += Buffer; 
 
Buffer = (-3.5297659265)*in[29]; 
out[6] += Buffer; 
 
Buffer = (-0.5473136903)*in[29]; 
out[7] += Buffer; 
 
Buffer = (-0.0188170046)*in[29]; 
out[8] += Buffer; 
 
Buffer = (-0.3901459189)*in[29]; 
out[9] += Buffer; 
 
Buffer = (-1.4719104438)*in[29]; 
out[10] += Buffer; 
 
Buffer = (-0.4110279866)*in[29]; 
out[11] += Buffer; 
 
Buffer = (4.1053953962)*in[29]; 
out[12] += Buffer; 
 
Buffer = (-0.3551221352)*in[29]; 
out[13] += Buffer; 
 
Buffer = (2.1653719455)*in[29]; 
out[14] += Buffer; 
 
Buffer = (-0.3003740274)*in[30]; 
out[0] += Buffer; 
 
Buffer = (-0.1762864376)*in[30]; 
out[1] += Buffer; 
 
Buffer = (0.0219980970)*in[30]; 
out[2] += Buffer; 
 
Buffer = (0.5238229173)*in[30]; 
out[3] += Buffer; 
 
Buffer = (-0.5473136903)*in[30]; 
out[4] += Buffer; 
 
Buffer = (-0.0188170046)*in[30]; 
out[5] += Buffer; 
 
Buffer = (0.1135037218)*in[30]; 
out[6] += Buffer; 
 
Buffer = (-3.5297659265)*in[30]; 
out[7] += Buffer; 
 
Buffer = (-0.3901459189)*in[30]; 
out[8] += Buffer; 
 
Buffer = (0.2706714932)*in[30]; 
out[9] += Buffer; 
 
Buffer = (-0.4110279866)*in[30]; 
out[10] += Buffer; 
 
Buffer = (4.1053953962)*in[30]; 
out[11] += Buffer; 
 
Buffer = (-1.4719104438)*in[30]; 
out[12] += Buffer; 
 
Buffer = (-0.3551221352)*in[30]; 
out[13] += Buffer; 
 
Buffer = (2.1653719455)*in[30]; 
out[14] += Buffer; 
 
Buffer = (-0.3003740274)*in[31]; 
out[0] += Buffer; 
 
Buffer = (0.0219980970)*in[31]; 
out[1] += Buffer; 
 
Buffer = (-0.1762864376)*in[31]; 
out[2] += Buffer; 
 
Buffer = (0.5238229173)*in[31]; 
out[3] += Buffer; 
 
Buffer = (0.1135037218)*in[31]; 
out[4] += Buffer; 
 
Buffer = (-0.0188170046)*in[31]; 
out[5] += Buffer; 
 
Buffer = (-0.5473136903)*in[31]; 
out[6] += Buffer; 
 
Buffer = (-3.5297659265)*in[31]; 
out[7] += Buffer; 
 
Buffer = (0.2706714932)*in[31]; 
out[8] += Buffer; 
 
Buffer = (-0.3901459189)*in[31]; 
out[9] += Buffer; 
 
Buffer = (-0.4110279866)*in[31]; 
out[10] += Buffer; 
 
Buffer = (-1.4719104438)*in[31]; 
out[11] += Buffer; 
 
Buffer = (4.1053953962)*in[31]; 
out[12] += Buffer; 
 
Buffer = (-0.3551221352)*in[31]; 
out[13] += Buffer; 
 
Buffer = (2.1653719455)*in[31]; 
out[14] += Buffer; 
 
Buffer = (-0.8241969448)*in[32]; 
out[0] += Buffer; 
 
Buffer = (0.8241969448)*in[32]; 
out[3] += Buffer; 
 
Buffer = (-0.1571677714)*in[32]; 
out[4] += Buffer; 
 
Buffer = (-0.1571677714)*in[32]; 
out[6] += Buffer; 
 
Buffer = (0.1571677714)*in[32]; 
out[8] += Buffer; 
 
Buffer = (0.1571677714)*in[32]; 
out[9] += Buffer; 
 
Buffer = (-0.0559058514)*in[32]; 
out[10] += Buffer; 
 
Buffer = (0.0559058514)*in[32]; 
out[13] += Buffer; 
 
Buffer = (0.1762864376)*in[33]; 
out[0] += Buffer; 
 
Buffer = (0.3003740274)*in[33]; 
out[1] += Buffer; 
 
Buffer = (-0.5238229173)*in[33]; 
out[2] += Buffer; 
 
Buffer = (-0.0219980970)*in[33]; 
out[3] += Buffer; 
 
Buffer = (0.5473136903)*in[33]; 
out[4] += Buffer; 
 
Buffer = (3.5297659265)*in[33]; 
out[5] += Buffer; 
 
Buffer = (0.3901459189)*in[33]; 
out[6] += Buffer; 
 
Buffer = (0.0188170046)*in[33]; 
out[7] += Buffer; 
 
Buffer = (-0.1135037218)*in[33]; 
out[8] += Buffer; 
 
Buffer = (-0.2706714932)*in[33]; 
out[9] += Buffer; 
 
Buffer = (-4.1053953962)*in[33]; 
out[10] += Buffer; 
 
Buffer = (0.4110279866)*in[33]; 
out[11] += Buffer; 
 
Buffer = (0.3551221352)*in[33]; 
out[12] += Buffer; 
 
Buffer = (1.4719104438)*in[33]; 
out[13] += Buffer; 
 
Buffer = (-2.1653719455)*in[33]; 
out[14] += Buffer; 
 
Buffer = (0.1762864376)*in[34]; 
out[0] += Buffer; 
 
Buffer = (-0.5238229173)*in[34]; 
out[1] += Buffer; 
 
Buffer = (0.3003740274)*in[34]; 
out[2] += Buffer; 
 
Buffer = (-0.0219980970)*in[34]; 
out[3] += Buffer; 
 
Buffer = (0.3901459189)*in[34]; 
out[4] += Buffer; 
 
Buffer = (3.5297659265)*in[34]; 
out[5] += Buffer; 
 
Buffer = (0.5473136903)*in[34]; 
out[6] += Buffer; 
 
Buffer = (0.0188170046)*in[34]; 
out[7] += Buffer; 
 
Buffer = (-0.2706714932)*in[34]; 
out[8] += Buffer; 
 
Buffer = (-0.1135037218)*in[34]; 
out[9] += Buffer; 
 
Buffer = (-4.1053953962)*in[34]; 
out[10] += Buffer; 
 
Buffer = (0.3551221352)*in[34]; 
out[11] += Buffer; 
 
Buffer = (0.4110279866)*in[34]; 
out[12] += Buffer; 
 
Buffer = (1.4719104438)*in[34]; 
out[13] += Buffer; 
 
Buffer = (-2.1653719455)*in[34]; 
out[14] += Buffer; 
 
Buffer = (0.1982845346)*in[35]; 
out[0] += Buffer; 
 
Buffer = (-0.1982845346)*in[35]; 
out[3] += Buffer; 
 
Buffer = (0.6608174121)*in[35]; 
out[4] += Buffer; 
 
Buffer = (0.6608174121)*in[35]; 
out[6] += Buffer; 
 
Buffer = (-0.6608174121)*in[35]; 
out[8] += Buffer; 
 
Buffer = (-0.6608174121)*in[35]; 
out[9] += Buffer; 
 
Buffer = (-5.5773058400)*in[35]; 
out[10] += Buffer; 
 
Buffer = (5.5773058400)*in[35]; 
out[13] += Buffer; 
 
Buffer = (0.1982845346)*in[36]; 
out[0] += Buffer; 
 
Buffer = (-0.1982845346)*in[36]; 
out[1] += Buffer; 
 
Buffer = (-0.6608174121)*in[36]; 
out[5] += Buffer; 
 
Buffer = (0.6608174121)*in[36]; 
out[6] += Buffer; 
 
Buffer = (0.6608174121)*in[36]; 
out[7] += Buffer; 
 
Buffer = (-0.6608174121)*in[36]; 
out[8] += Buffer; 
 
Buffer = (-5.5773058400)*in[36]; 
out[12] += Buffer; 
 
Buffer = (5.5773058400)*in[36]; 
out[13] += Buffer; 
 
Buffer = (0.1762864376)*in[37]; 
out[0] += Buffer; 
 
Buffer = (-0.0219980970)*in[37]; 
out[1] += Buffer; 
 
Buffer = (0.3003740274)*in[37]; 
out[2] += Buffer; 
 
Buffer = (-0.5238229173)*in[37]; 
out[3] += Buffer; 
 
Buffer = (0.0188170046)*in[37]; 
out[4] += Buffer; 
 
Buffer = (-0.1135037218)*in[37]; 
out[5] += Buffer; 
 
Buffer = (0.5473136903)*in[37]; 
out[6] += Buffer; 
 
Buffer = (0.3901459189)*in[37]; 
out[7] += Buffer; 
 
Buffer = (-0.2706714932)*in[37]; 
out[8] += Buffer; 
 
Buffer = (3.5297659265)*in[37]; 
out[9] += Buffer; 
 
Buffer = (0.4110279866)*in[37]; 
out[10] += Buffer; 
 
Buffer = (0.3551221352)*in[37]; 
out[11] += Buffer; 
 
Buffer = (-4.1053953962)*in[37]; 
out[12] += Buffer; 
 
Buffer = (1.4719104438)*in[37]; 
out[13] += Buffer; 
 
Buffer = (-2.1653719455)*in[37]; 
out[14] += Buffer; 
 
Buffer = (0.1762864376)*in[38]; 
out[0] += Buffer; 
 
Buffer = (-0.0219980970)*in[38]; 
out[1] += Buffer; 
 
Buffer = (-0.5238229173)*in[38]; 
out[2] += Buffer; 
 
Buffer = (0.3003740274)*in[38]; 
out[3] += Buffer; 
 
Buffer = (0.0188170046)*in[38]; 
out[4] += Buffer; 
 
Buffer = (-0.2706714932)*in[38]; 
out[5] += Buffer; 
 
Buffer = (0.3901459189)*in[38]; 
out[6] += Buffer; 
 
Buffer = (0.5473136903)*in[38]; 
out[7] += Buffer; 
 
Buffer = (-0.1135037218)*in[38]; 
out[8] += Buffer; 
 
Buffer = (3.5297659265)*in[38]; 
out[9] += Buffer; 
 
Buffer = (0.3551221352)*in[38]; 
out[10] += Buffer; 
 
Buffer = (0.4110279866)*in[38]; 
out[11] += Buffer; 
 
Buffer = (-4.1053953962)*in[38]; 
out[12] += Buffer; 
 
Buffer = (1.4719104438)*in[38]; 
out[13] += Buffer; 
 
Buffer = (-2.1653719455)*in[38]; 
out[14] += Buffer; 
 
Buffer = (0.1762864376)*in[39]; 
out[0] += Buffer; 
 
Buffer = (0.3003740274)*in[39]; 
out[1] += Buffer; 
 
Buffer = (-0.0219980970)*in[39]; 
out[2] += Buffer; 
 
Buffer = (-0.5238229173)*in[39]; 
out[3] += Buffer; 
 
Buffer = (0.5473136903)*in[39]; 
out[4] += Buffer; 
 
Buffer = (-0.1135037218)*in[39]; 
out[5] += Buffer; 
 
Buffer = (0.0188170046)*in[39]; 
out[6] += Buffer; 
 
Buffer = (0.3901459189)*in[39]; 
out[7] += Buffer; 
 
Buffer = (3.5297659265)*in[39]; 
out[8] += Buffer; 
 
Buffer = (-0.2706714932)*in[39]; 
out[9] += Buffer; 
 
Buffer = (0.4110279866)*in[39]; 
out[10] += Buffer; 
 
Buffer = (-4.1053953962)*in[39]; 
out[11] += Buffer; 
 
Buffer = (0.3551221352)*in[39]; 
out[12] += Buffer; 
 
Buffer = (1.4719104438)*in[39]; 
out[13] += Buffer; 
 
Buffer = (-2.1653719455)*in[39]; 
out[14] += Buffer; 
 
Buffer = (0.1982845346)*in[40]; 
out[0] += Buffer; 
 
Buffer = (-0.1982845346)*in[40]; 
out[2] += Buffer; 
 
Buffer = (0.6608174121)*in[40]; 
out[4] += Buffer; 
 
Buffer = (-0.6608174121)*in[40]; 
out[5] += Buffer; 
 
Buffer = (0.6608174121)*in[40]; 
out[7] += Buffer; 
 
Buffer = (-0.6608174121)*in[40]; 
out[9] += Buffer; 
 
Buffer = (-5.5773058400)*in[40]; 
out[11] += Buffer; 
 
Buffer = (5.5773058400)*in[40]; 
out[13] += Buffer; 
 
Buffer = (0.1762864376)*in[41]; 
out[0] += Buffer; 
 
Buffer = (-0.5238229173)*in[41]; 
out[1] += Buffer; 
 
Buffer = (-0.0219980970)*in[41]; 
out[2] += Buffer; 
 
Buffer = (0.3003740274)*in[41]; 
out[3] += Buffer; 
 
Buffer = (0.3901459189)*in[41]; 
out[4] += Buffer; 
 
Buffer = (-0.2706714932)*in[41]; 
out[5] += Buffer; 
 
Buffer = (0.0188170046)*in[41]; 
out[6] += Buffer; 
 
Buffer = (0.5473136903)*in[41]; 
out[7] += Buffer; 
 
Buffer = (3.5297659265)*in[41]; 
out[8] += Buffer; 
 
Buffer = (-0.1135037218)*in[41]; 
out[9] += Buffer; 
 
Buffer = (0.3551221352)*in[41]; 
out[10] += Buffer; 
 
Buffer = (-4.1053953962)*in[41]; 
out[11] += Buffer; 
 
Buffer = (0.4110279866)*in[41]; 
out[12] += Buffer; 
 
Buffer = (1.4719104438)*in[41]; 
out[13] += Buffer; 
 
Buffer = (-2.1653719455)*in[41]; 
out[14] += Buffer; 
 



            };
            
            inline void GradientInterpolation(const Real * in, Real *out)
            {
                const Real c1 = TetriX::Tet_14::_c1;
                const Real c2 = TetriX::Tet_14::_c2;
                const Real d = TetriX::Tet_14::_d;

                const Real c1_2 = c1*c1;
                const Real c2_2 = c2*c2;
                const Real d_2 = d*d;

                const Real c1_3 = c1*c1*c1;
                const Real c2_3 = c2*c2*c2;
                const Real d_3 = d*d*d;

                Real Buffer;

                for (Index i=0;i<(3*14);++i) 
                {

                    out[i]=0.0;
                }

Buffer = (-1.5846000073)*in[0]; 
out[0] += Buffer; 
 
Buffer = (-1.5846000073)*in[0]; 
out[1] += Buffer; 
 
Buffer = (-1.5846000073)*in[0]; 
out[2] += Buffer; 
 
Buffer = (0.2748830499)*in[0]; 
out[3] += Buffer; 
 
Buffer = (0.4282516877)*in[0]; 
out[4] += Buffer; 
 
Buffer = (0.4282516877)*in[0]; 
out[5] += Buffer; 
 
Buffer = (0.4282516877)*in[0]; 
out[6] += Buffer; 
 
Buffer = (0.2748830499)*in[0]; 
out[7] += Buffer; 
 
Buffer = (0.4282516877)*in[0]; 
out[8] += Buffer; 
 
Buffer = (0.4282516877)*in[0]; 
out[9] += Buffer; 
 
Buffer = (0.4282516877)*in[0]; 
out[10] += Buffer; 
 
Buffer = (0.2748830499)*in[0]; 
out[11] += Buffer; 
 
Buffer = (0.0805491062)*in[0]; 
out[12] += Buffer; 
 
Buffer = (0.0805491062)*in[0]; 
out[13] += Buffer; 
 
Buffer = (0.0805491062)*in[0]; 
out[14] += Buffer; 
 
Buffer = (-0.1733620808)*in[0]; 
out[15] += Buffer; 
 
Buffer = (-0.3063509417)*in[0]; 
out[16] += Buffer; 
 
Buffer = (-0.3063509417)*in[0]; 
out[17] += Buffer; 
 
Buffer = (-0.3063509417)*in[0]; 
out[18] += Buffer; 
 
Buffer = (-0.1733620808)*in[0]; 
out[19] += Buffer; 
 
Buffer = (-0.3063509417)*in[0]; 
out[20] += Buffer; 
 
Buffer = (-0.3063509417)*in[0]; 
out[21] += Buffer; 
 
Buffer = (-0.3063509417)*in[0]; 
out[22] += Buffer; 
 
Buffer = (-0.1733620808)*in[0]; 
out[23] += Buffer; 
 
Buffer = (-0.8241969448)*in[0]; 
out[24] += Buffer; 
 
Buffer = (-0.3003740274)*in[0]; 
out[25] += Buffer; 
 
Buffer = (-0.3003740274)*in[0]; 
out[26] += Buffer; 
 
Buffer = (-0.3003740274)*in[0]; 
out[27] += Buffer; 
 
Buffer = (-0.8241969448)*in[0]; 
out[28] += Buffer; 
 
Buffer = (-0.3003740274)*in[0]; 
out[29] += Buffer; 
 
Buffer = (-0.3003740274)*in[0]; 
out[30] += Buffer; 
 
Buffer = (-0.3003740274)*in[0]; 
out[31] += Buffer; 
 
Buffer = (-0.8241969448)*in[0]; 
out[32] += Buffer; 
 
Buffer = (0.1762864376)*in[0]; 
out[33] += Buffer; 
 
Buffer = (0.1762864376)*in[0]; 
out[34] += Buffer; 
 
Buffer = (0.1982845346)*in[0]; 
out[35] += Buffer; 
 
Buffer = (0.1982845346)*in[0]; 
out[36] += Buffer; 
 
Buffer = (0.1762864376)*in[0]; 
out[37] += Buffer; 
 
Buffer = (0.1762864376)*in[0]; 
out[38] += Buffer; 
 
Buffer = (0.1762864376)*in[0]; 
out[39] += Buffer; 
 
Buffer = (0.1982845346)*in[0]; 
out[40] += Buffer; 
 
Buffer = (0.1762864376)*in[0]; 
out[41] += Buffer; 
 
Buffer = (-0.2748830499)*in[1]; 
out[0] += Buffer; 
 
Buffer = (0.1533686378)*in[1]; 
out[1] += Buffer; 
 
Buffer = (0.1533686378)*in[1]; 
out[2] += Buffer; 
 
Buffer = (1.5846000073)*in[1]; 
out[3] += Buffer; 
 
Buffer = (-0.4282516877)*in[1]; 
out[6] += Buffer; 
 
Buffer = (-0.1533686378)*in[1]; 
out[7] += Buffer; 
 
Buffer = (-0.4282516877)*in[1]; 
out[9] += Buffer; 
 
Buffer = (-0.1533686378)*in[1]; 
out[11] += Buffer; 
 
Buffer = (0.1733620808)*in[1]; 
out[12] += Buffer; 
 
Buffer = (-0.1329888608)*in[1]; 
out[13] += Buffer; 
 
Buffer = (-0.1329888608)*in[1]; 
out[14] += Buffer; 
 
Buffer = (-0.0805491062)*in[1]; 
out[15] += Buffer; 
 
Buffer = (0.3063509417)*in[1]; 
out[18] += Buffer; 
 
Buffer = (0.1329888608)*in[1]; 
out[19] += Buffer; 
 
Buffer = (0.3063509417)*in[1]; 
out[21] += Buffer; 
 
Buffer = (0.1329888608)*in[1]; 
out[23] += Buffer; 
 
Buffer = (0.8241969448)*in[1]; 
out[24] += Buffer; 
 
Buffer = (0.5238229173)*in[1]; 
out[25] += Buffer; 
 
Buffer = (0.5238229173)*in[1]; 
out[26] += Buffer; 
 
Buffer = (-0.1762864376)*in[1]; 
out[27] += Buffer; 
 
Buffer = (0.0219980970)*in[1]; 
out[29] += Buffer; 
 
Buffer = (-0.1762864376)*in[1]; 
out[30] += Buffer; 
 
Buffer = (0.0219980970)*in[1]; 
out[31] += Buffer; 
 
Buffer = (0.3003740274)*in[1]; 
out[33] += Buffer; 
 
Buffer = (-0.5238229173)*in[1]; 
out[34] += Buffer; 
 
Buffer = (-0.1982845346)*in[1]; 
out[36] += Buffer; 
 
Buffer = (-0.0219980970)*in[1]; 
out[37] += Buffer; 
 
Buffer = (-0.0219980970)*in[1]; 
out[38] += Buffer; 
 
Buffer = (0.3003740274)*in[1]; 
out[39] += Buffer; 
 
Buffer = (-0.5238229173)*in[1]; 
out[41] += Buffer; 
 
Buffer = (0.1533686378)*in[2]; 
out[0] += Buffer; 
 
Buffer = (-0.2748830499)*in[2]; 
out[1] += Buffer; 
 
Buffer = (0.1533686378)*in[2]; 
out[2] += Buffer; 
 
Buffer = (-0.1533686378)*in[2]; 
out[3] += Buffer; 
 
Buffer = (-0.4282516877)*in[2]; 
out[4] += Buffer; 
 
Buffer = (1.5846000073)*in[2]; 
out[7] += Buffer; 
 
Buffer = (-0.4282516877)*in[2]; 
out[10] += Buffer; 
 
Buffer = (-0.1533686378)*in[2]; 
out[11] += Buffer; 
 
Buffer = (-0.1329888608)*in[2]; 
out[12] += Buffer; 
 
Buffer = (0.1733620808)*in[2]; 
out[13] += Buffer; 
 
Buffer = (-0.1329888608)*in[2]; 
out[14] += Buffer; 
 
Buffer = (0.1329888608)*in[2]; 
out[15] += Buffer; 
 
Buffer = (0.3063509417)*in[2]; 
out[16] += Buffer; 
 
Buffer = (-0.0805491062)*in[2]; 
out[19] += Buffer; 
 
Buffer = (0.3063509417)*in[2]; 
out[22] += Buffer; 
 
Buffer = (0.1329888608)*in[2]; 
out[23] += Buffer; 
 
Buffer = (-0.1762864376)*in[2]; 
out[25] += Buffer; 
 
Buffer = (0.0219980970)*in[2]; 
out[26] += Buffer; 
 
Buffer = (0.5238229173)*in[2]; 
out[27] += Buffer; 
 
Buffer = (0.8241969448)*in[2]; 
out[28] += Buffer; 
 
Buffer = (0.5238229173)*in[2]; 
out[29] += Buffer; 
 
Buffer = (0.0219980970)*in[2]; 
out[30] += Buffer; 
 
Buffer = (-0.1762864376)*in[2]; 
out[31] += Buffer; 
 
Buffer = (-0.5238229173)*in[2]; 
out[33] += Buffer; 
 
Buffer = (0.3003740274)*in[2]; 
out[34] += Buffer; 
 
Buffer = (0.3003740274)*in[2]; 
out[37] += Buffer; 
 
Buffer = (-0.5238229173)*in[2]; 
out[38] += Buffer; 
 
Buffer = (-0.0219980970)*in[2]; 
out[39] += Buffer; 
 
Buffer = (-0.1982845346)*in[2]; 
out[40] += Buffer; 
 
Buffer = (-0.0219980970)*in[2]; 
out[41] += Buffer; 
 
Buffer = (0.1533686378)*in[3]; 
out[0] += Buffer; 
 
Buffer = (0.1533686378)*in[3]; 
out[1] += Buffer; 
 
Buffer = (-0.2748830499)*in[3]; 
out[2] += Buffer; 
 
Buffer = (-0.1533686378)*in[3]; 
out[3] += Buffer; 
 
Buffer = (-0.4282516877)*in[3]; 
out[5] += Buffer; 
 
Buffer = (-0.1533686378)*in[3]; 
out[7] += Buffer; 
 
Buffer = (-0.4282516877)*in[3]; 
out[8] += Buffer; 
 
Buffer = (1.5846000073)*in[3]; 
out[11] += Buffer; 
 
Buffer = (-0.1329888608)*in[3]; 
out[12] += Buffer; 
 
Buffer = (-0.1329888608)*in[3]; 
out[13] += Buffer; 
 
Buffer = (0.1733620808)*in[3]; 
out[14] += Buffer; 
 
Buffer = (0.1329888608)*in[3]; 
out[15] += Buffer; 
 
Buffer = (0.3063509417)*in[3]; 
out[17] += Buffer; 
 
Buffer = (0.1329888608)*in[3]; 
out[19] += Buffer; 
 
Buffer = (0.3063509417)*in[3]; 
out[20] += Buffer; 
 
Buffer = (-0.0805491062)*in[3]; 
out[23] += Buffer; 
 
Buffer = (0.0219980970)*in[3]; 
out[25] += Buffer; 
 
Buffer = (-0.1762864376)*in[3]; 
out[26] += Buffer; 
 
Buffer = (0.0219980970)*in[3]; 
out[27] += Buffer; 
 
Buffer = (-0.1762864376)*in[3]; 
out[29] += Buffer; 
 
Buffer = (0.5238229173)*in[3]; 
out[30] += Buffer; 
 
Buffer = (0.5238229173)*in[3]; 
out[31] += Buffer; 
 
Buffer = (0.8241969448)*in[3]; 
out[32] += Buffer; 
 
Buffer = (-0.0219980970)*in[3]; 
out[33] += Buffer; 
 
Buffer = (-0.0219980970)*in[3]; 
out[34] += Buffer; 
 
Buffer = (-0.1982845346)*in[3]; 
out[35] += Buffer; 
 
Buffer = (-0.5238229173)*in[3]; 
out[37] += Buffer; 
 
Buffer = (0.3003740274)*in[3]; 
out[38] += Buffer; 
 
Buffer = (-0.5238229173)*in[3]; 
out[39] += Buffer; 
 
Buffer = (0.3003740274)*in[3]; 
out[41] += Buffer; 
 
Buffer = (1.2892868929)*in[4]; 
out[0] += Buffer; 
 
Buffer = (-0.7946608590)*in[4]; 
out[1] += Buffer; 
 
Buffer = (-0.7946608590)*in[4]; 
out[2] += Buffer; 
 
Buffer = (-1.2892868929)*in[4]; 
out[3] += Buffer; 
 
Buffer = (-2.0839477519)*in[4]; 
out[4] += Buffer; 
 
Buffer = (-2.0839477519)*in[4]; 
out[5] += Buffer; 
 
Buffer = (0.2591746961)*in[4]; 
out[7] += Buffer; 
 
Buffer = (0.4322882431)*in[4]; 
out[8] += Buffer; 
 
Buffer = (0.4322882431)*in[4]; 
out[10] += Buffer; 
 
Buffer = (0.2591746961)*in[4]; 
out[11] += Buffer; 
 
Buffer = (0.0897361785)*in[4]; 
out[12] += Buffer; 
 
Buffer = (0.0715962681)*in[4]; 
out[13] += Buffer; 
 
Buffer = (0.0715962681)*in[4]; 
out[14] += Buffer; 
 
Buffer = (-0.0897361785)*in[4]; 
out[15] += Buffer; 
 
Buffer = (-0.0181399103)*in[4]; 
out[16] += Buffer; 
 
Buffer = (-0.0181399103)*in[4]; 
out[17] += Buffer; 
 
Buffer = (-0.2390823071)*in[4]; 
out[19] += Buffer; 
 
Buffer = (-0.9923146186)*in[4]; 
out[20] += Buffer; 
 
Buffer = (-0.9923146186)*in[4]; 
out[22] += Buffer; 
 
Buffer = (-0.2390823071)*in[4]; 
out[23] += Buffer; 
 
Buffer = (-3.5297659265)*in[4]; 
out[25] += Buffer; 
 
Buffer = (-3.5297659265)*in[4]; 
out[26] += Buffer; 
 
Buffer = (-0.5473136903)*in[4]; 
out[27] += Buffer; 
 
Buffer = (-0.1571677714)*in[4]; 
out[28] += Buffer; 
 
Buffer = (0.1135037218)*in[4]; 
out[29] += Buffer; 
 
Buffer = (-0.5473136903)*in[4]; 
out[30] += Buffer; 
 
Buffer = (0.1135037218)*in[4]; 
out[31] += Buffer; 
 
Buffer = (-0.1571677714)*in[4]; 
out[32] += Buffer; 
 
Buffer = (0.5473136903)*in[4]; 
out[33] += Buffer; 
 
Buffer = (0.3901459189)*in[4]; 
out[34] += Buffer; 
 
Buffer = (0.6608174121)*in[4]; 
out[35] += Buffer; 
 
Buffer = (0.0188170046)*in[4]; 
out[37] += Buffer; 
 
Buffer = (0.0188170046)*in[4]; 
out[38] += Buffer; 
 
Buffer = (0.5473136903)*in[4]; 
out[39] += Buffer; 
 
Buffer = (0.6608174121)*in[4]; 
out[40] += Buffer; 
 
Buffer = (0.3901459189)*in[4]; 
out[41] += Buffer; 
 
Buffer = (-0.2591746961)*in[5]; 
out[0] += Buffer; 
 
Buffer = (-0.2591746961)*in[5]; 
out[1] += Buffer; 
 
Buffer = (0.1731135469)*in[5]; 
out[2] += Buffer; 
 
Buffer = (0.7946608590)*in[5]; 
out[3] += Buffer; 
 
Buffer = (2.0839477519)*in[5]; 
out[4] += Buffer; 
 
Buffer = (2.0839477519)*in[5]; 
out[6] += Buffer; 
 
Buffer = (0.7946608590)*in[5]; 
out[7] += Buffer; 
 
Buffer = (-0.4322882431)*in[5]; 
out[9] += Buffer; 
 
Buffer = (-0.4322882431)*in[5]; 
out[10] += Buffer; 
 
Buffer = (-0.1731135469)*in[5]; 
out[11] += Buffer; 
 
Buffer = (0.2390823071)*in[5]; 
out[12] += Buffer; 
 
Buffer = (0.2390823071)*in[5]; 
out[13] += Buffer; 
 
Buffer = (-0.7532323114)*in[5]; 
out[14] += Buffer; 
 
Buffer = (-0.0715962681)*in[5]; 
out[15] += Buffer; 
 
Buffer = (0.0181399103)*in[5]; 
out[16] += Buffer; 
 
Buffer = (0.0181399103)*in[5]; 
out[18] += Buffer; 
 
Buffer = (-0.0715962681)*in[5]; 
out[19] += Buffer; 
 
Buffer = (0.9923146186)*in[5]; 
out[21] += Buffer; 
 
Buffer = (0.9923146186)*in[5]; 
out[22] += Buffer; 
 
Buffer = (0.7532323114)*in[5]; 
out[23] += Buffer; 
 
Buffer = (0.1571677714)*in[5]; 
out[24] += Buffer; 
 
Buffer = (-0.3901459189)*in[5]; 
out[25] += Buffer; 
 
Buffer = (0.2706714932)*in[5]; 
out[26] += Buffer; 
 
Buffer = (-0.3901459189)*in[5]; 
out[27] += Buffer; 
 
Buffer = (0.1571677714)*in[5]; 
out[28] += Buffer; 
 
Buffer = (0.2706714932)*in[5]; 
out[29] += Buffer; 
 
Buffer = (-0.0188170046)*in[5]; 
out[30] += Buffer; 
 
Buffer = (-0.0188170046)*in[5]; 
out[31] += Buffer; 
 
Buffer = (3.5297659265)*in[5]; 
out[33] += Buffer; 
 
Buffer = (3.5297659265)*in[5]; 
out[34] += Buffer; 
 
Buffer = (-0.6608174121)*in[5]; 
out[36] += Buffer; 
 
Buffer = (-0.1135037218)*in[5]; 
out[37] += Buffer; 
 
Buffer = (-0.2706714932)*in[5]; 
out[38] += Buffer; 
 
Buffer = (-0.1135037218)*in[5]; 
out[39] += Buffer; 
 
Buffer = (-0.6608174121)*in[5]; 
out[40] += Buffer; 
 
Buffer = (-0.2706714932)*in[5]; 
out[41] += Buffer; 
 
Buffer = (-0.7946608590)*in[6]; 
out[0] += Buffer; 
 
Buffer = (1.2892868929)*in[6]; 
out[1] += Buffer; 
 
Buffer = (-0.7946608590)*in[6]; 
out[2] += Buffer; 
 
Buffer = (0.2591746961)*in[6]; 
out[3] += Buffer; 
 
Buffer = (0.4322882431)*in[6]; 
out[5] += Buffer; 
 
Buffer = (-2.0839477519)*in[6]; 
out[6] += Buffer; 
 
Buffer = (-1.2892868929)*in[6]; 
out[7] += Buffer; 
 
Buffer = (-2.0839477519)*in[6]; 
out[8] += Buffer; 
 
Buffer = (0.4322882431)*in[6]; 
out[9] += Buffer; 
 
Buffer = (0.2591746961)*in[6]; 
out[11] += Buffer; 
 
Buffer = (0.0715962681)*in[6]; 
out[12] += Buffer; 
 
Buffer = (0.0897361785)*in[6]; 
out[13] += Buffer; 
 
Buffer = (0.0715962681)*in[6]; 
out[14] += Buffer; 
 
Buffer = (-0.2390823071)*in[6]; 
out[15] += Buffer; 
 
Buffer = (-0.9923146186)*in[6]; 
out[17] += Buffer; 
 
Buffer = (-0.0181399103)*in[6]; 
out[18] += Buffer; 
 
Buffer = (-0.0897361785)*in[6]; 
out[19] += Buffer; 
 
Buffer = (-0.0181399103)*in[6]; 
out[20] += Buffer; 
 
Buffer = (-0.9923146186)*in[6]; 
out[21] += Buffer; 
 
Buffer = (-0.2390823071)*in[6]; 
out[23] += Buffer; 
 
Buffer = (-0.1571677714)*in[6]; 
out[24] += Buffer; 
 
Buffer = (-0.5473136903)*in[6]; 
out[25] += Buffer; 
 
Buffer = (0.1135037218)*in[6]; 
out[26] += Buffer; 
 
Buffer = (-3.5297659265)*in[6]; 
out[27] += Buffer; 
 
Buffer = (-3.5297659265)*in[6]; 
out[29] += Buffer; 
 
Buffer = (0.1135037218)*in[6]; 
out[30] += Buffer; 
 
Buffer = (-0.5473136903)*in[6]; 
out[31] += Buffer; 
 
Buffer = (-0.1571677714)*in[6]; 
out[32] += Buffer; 
 
Buffer = (0.3901459189)*in[6]; 
out[33] += Buffer; 
 
Buffer = (0.5473136903)*in[6]; 
out[34] += Buffer; 
 
Buffer = (0.6608174121)*in[6]; 
out[35] += Buffer; 
 
Buffer = (0.6608174121)*in[6]; 
out[36] += Buffer; 
 
Buffer = (0.5473136903)*in[6]; 
out[37] += Buffer; 
 
Buffer = (0.3901459189)*in[6]; 
out[38] += Buffer; 
 
Buffer = (0.0188170046)*in[6]; 
out[39] += Buffer; 
 
Buffer = (0.0188170046)*in[6]; 
out[41] += Buffer; 
 
Buffer = (-0.7946608590)*in[7]; 
out[0] += Buffer; 
 
Buffer = (-0.7946608590)*in[7]; 
out[1] += Buffer; 
 
Buffer = (1.2892868929)*in[7]; 
out[2] += Buffer; 
 
Buffer = (0.2591746961)*in[7]; 
out[3] += Buffer; 
 
Buffer = (0.4322882431)*in[7]; 
out[4] += Buffer; 
 
Buffer = (0.4322882431)*in[7]; 
out[6] += Buffer; 
 
Buffer = (0.2591746961)*in[7]; 
out[7] += Buffer; 
 
Buffer = (-2.0839477519)*in[7]; 
out[9] += Buffer; 
 
Buffer = (-2.0839477519)*in[7]; 
out[10] += Buffer; 
 
Buffer = (-1.2892868929)*in[7]; 
out[11] += Buffer; 
 
Buffer = (0.0715962681)*in[7]; 
out[12] += Buffer; 
 
Buffer = (0.0715962681)*in[7]; 
out[13] += Buffer; 
 
Buffer = (0.0897361785)*in[7]; 
out[14] += Buffer; 
 
Buffer = (-0.2390823071)*in[7]; 
out[15] += Buffer; 
 
Buffer = (-0.9923146186)*in[7]; 
out[16] += Buffer; 
 
Buffer = (-0.9923146186)*in[7]; 
out[18] += Buffer; 
 
Buffer = (-0.2390823071)*in[7]; 
out[19] += Buffer; 
 
Buffer = (-0.0181399103)*in[7]; 
out[21] += Buffer; 
 
Buffer = (-0.0181399103)*in[7]; 
out[22] += Buffer; 
 
Buffer = (-0.0897361785)*in[7]; 
out[23] += Buffer; 
 
Buffer = (-0.1571677714)*in[7]; 
out[24] += Buffer; 
 
Buffer = (0.1135037218)*in[7]; 
out[25] += Buffer; 
 
Buffer = (-0.5473136903)*in[7]; 
out[26] += Buffer; 
 
Buffer = (0.1135037218)*in[7]; 
out[27] += Buffer; 
 
Buffer = (-0.1571677714)*in[7]; 
out[28] += Buffer; 
 
Buffer = (-0.5473136903)*in[7]; 
out[29] += Buffer; 
 
Buffer = (-3.5297659265)*in[7]; 
out[30] += Buffer; 
 
Buffer = (-3.5297659265)*in[7]; 
out[31] += Buffer; 
 
Buffer = (0.0188170046)*in[7]; 
out[33] += Buffer; 
 
Buffer = (0.0188170046)*in[7]; 
out[34] += Buffer; 
 
Buffer = (0.6608174121)*in[7]; 
out[36] += Buffer; 
 
Buffer = (0.3901459189)*in[7]; 
out[37] += Buffer; 
 
Buffer = (0.5473136903)*in[7]; 
out[38] += Buffer; 
 
Buffer = (0.3901459189)*in[7]; 
out[39] += Buffer; 
 
Buffer = (0.6608174121)*in[7]; 
out[40] += Buffer; 
 
Buffer = (0.5473136903)*in[7]; 
out[41] += Buffer; 
 
Buffer = (-0.2591746961)*in[8]; 
out[0] += Buffer; 
 
Buffer = (0.1731135469)*in[8]; 
out[1] += Buffer; 
 
Buffer = (-0.2591746961)*in[8]; 
out[2] += Buffer; 
 
Buffer = (0.7946608590)*in[8]; 
out[3] += Buffer; 
 
Buffer = (2.0839477519)*in[8]; 
out[5] += Buffer; 
 
Buffer = (-0.4322882431)*in[8]; 
out[6] += Buffer; 
 
Buffer = (-0.1731135469)*in[8]; 
out[7] += Buffer; 
 
Buffer = (-0.4322882431)*in[8]; 
out[8] += Buffer; 
 
Buffer = (2.0839477519)*in[8]; 
out[9] += Buffer; 
 
Buffer = (0.7946608590)*in[8]; 
out[11] += Buffer; 
 
Buffer = (0.2390823071)*in[8]; 
out[12] += Buffer; 
 
Buffer = (-0.7532323114)*in[8]; 
out[13] += Buffer; 
 
Buffer = (0.2390823071)*in[8]; 
out[14] += Buffer; 
 
Buffer = (-0.0715962681)*in[8]; 
out[15] += Buffer; 
 
Buffer = (0.0181399103)*in[8]; 
out[17] += Buffer; 
 
Buffer = (0.9923146186)*in[8]; 
out[18] += Buffer; 
 
Buffer = (0.7532323114)*in[8]; 
out[19] += Buffer; 
 
Buffer = (0.9923146186)*in[8]; 
out[20] += Buffer; 
 
Buffer = (0.0181399103)*in[8]; 
out[21] += Buffer; 
 
Buffer = (-0.0715962681)*in[8]; 
out[23] += Buffer; 
 
Buffer = (0.1571677714)*in[8]; 
out[24] += Buffer; 
 
Buffer = (0.2706714932)*in[8]; 
out[25] += Buffer; 
 
Buffer = (-0.3901459189)*in[8]; 
out[26] += Buffer; 
 
Buffer = (-0.0188170046)*in[8]; 
out[27] += Buffer; 
 
Buffer = (-0.0188170046)*in[8]; 
out[29] += Buffer; 
 
Buffer = (-0.3901459189)*in[8]; 
out[30] += Buffer; 
 
Buffer = (0.2706714932)*in[8]; 
out[31] += Buffer; 
 
Buffer = (0.1571677714)*in[8]; 
out[32] += Buffer; 
 
Buffer = (-0.1135037218)*in[8]; 
out[33] += Buffer; 
 
Buffer = (-0.2706714932)*in[8]; 
out[34] += Buffer; 
 
Buffer = (-0.6608174121)*in[8]; 
out[35] += Buffer; 
 
Buffer = (-0.6608174121)*in[8]; 
out[36] += Buffer; 
 
Buffer = (-0.2706714932)*in[8]; 
out[37] += Buffer; 
 
Buffer = (-0.1135037218)*in[8]; 
out[38] += Buffer; 
 
Buffer = (3.5297659265)*in[8]; 
out[39] += Buffer; 
 
Buffer = (3.5297659265)*in[8]; 
out[41] += Buffer; 
 
Buffer = (0.1731135469)*in[9]; 
out[0] += Buffer; 
 
Buffer = (-0.2591746961)*in[9]; 
out[1] += Buffer; 
 
Buffer = (-0.2591746961)*in[9]; 
out[2] += Buffer; 
 
Buffer = (-0.1731135469)*in[9]; 
out[3] += Buffer; 
 
Buffer = (-0.4322882431)*in[9]; 
out[4] += Buffer; 
 
Buffer = (-0.4322882431)*in[9]; 
out[5] += Buffer; 
 
Buffer = (0.7946608590)*in[9]; 
out[7] += Buffer; 
 
Buffer = (2.0839477519)*in[9]; 
out[8] += Buffer; 
 
Buffer = (2.0839477519)*in[9]; 
out[10] += Buffer; 
 
Buffer = (0.7946608590)*in[9]; 
out[11] += Buffer; 
 
Buffer = (-0.7532323114)*in[9]; 
out[12] += Buffer; 
 
Buffer = (0.2390823071)*in[9]; 
out[13] += Buffer; 
 
Buffer = (0.2390823071)*in[9]; 
out[14] += Buffer; 
 
Buffer = (0.7532323114)*in[9]; 
out[15] += Buffer; 
 
Buffer = (0.9923146186)*in[9]; 
out[16] += Buffer; 
 
Buffer = (0.9923146186)*in[9]; 
out[17] += Buffer; 
 
Buffer = (-0.0715962681)*in[9]; 
out[19] += Buffer; 
 
Buffer = (0.0181399103)*in[9]; 
out[20] += Buffer; 
 
Buffer = (0.0181399103)*in[9]; 
out[22] += Buffer; 
 
Buffer = (-0.0715962681)*in[9]; 
out[23] += Buffer; 
 
Buffer = (-0.0188170046)*in[9]; 
out[25] += Buffer; 
 
Buffer = (-0.0188170046)*in[9]; 
out[26] += Buffer; 
 
Buffer = (0.2706714932)*in[9]; 
out[27] += Buffer; 
 
Buffer = (0.1571677714)*in[9]; 
out[28] += Buffer; 
 
Buffer = (-0.3901459189)*in[9]; 
out[29] += Buffer; 
 
Buffer = (0.2706714932)*in[9]; 
out[30] += Buffer; 
 
Buffer = (-0.3901459189)*in[9]; 
out[31] += Buffer; 
 
Buffer = (0.1571677714)*in[9]; 
out[32] += Buffer; 
 
Buffer = (-0.2706714932)*in[9]; 
out[33] += Buffer; 
 
Buffer = (-0.1135037218)*in[9]; 
out[34] += Buffer; 
 
Buffer = (-0.6608174121)*in[9]; 
out[35] += Buffer; 
 
Buffer = (3.5297659265)*in[9]; 
out[37] += Buffer; 
 
Buffer = (3.5297659265)*in[9]; 
out[38] += Buffer; 
 
Buffer = (-0.2706714932)*in[9]; 
out[39] += Buffer; 
 
Buffer = (-0.6608174121)*in[9]; 
out[40] += Buffer; 
 
Buffer = (-0.1135037218)*in[9]; 
out[41] += Buffer; 
 
Buffer = (0.9908122592)*in[10]; 
out[0] += Buffer; 
 
Buffer = (0.9908122592)*in[10]; 
out[1] += Buffer; 
 
Buffer = (-0.8164535405)*in[10]; 
out[2] += Buffer; 
 
Buffer = (-0.9908122592)*in[10]; 
out[3] += Buffer; 
 
Buffer = (-1.8072657997)*in[10]; 
out[5] += Buffer; 
 
Buffer = (-0.9908122592)*in[10]; 
out[7] += Buffer; 
 
Buffer = (-1.8072657997)*in[10]; 
out[8] += Buffer; 
 
Buffer = (0.3520629014)*in[10]; 
out[11] += Buffer; 
 
Buffer = (0.4978729532)*in[10]; 
out[12] += Buffer; 
 
Buffer = (0.4978729532)*in[10]; 
out[13] += Buffer; 
 
Buffer = (-0.0673924284)*in[10]; 
out[14] += Buffer; 
 
Buffer = (-0.4978729532)*in[10]; 
out[15] += Buffer; 
 
Buffer = (-0.5652653816)*in[10]; 
out[17] += Buffer; 
 
Buffer = (-0.4978729532)*in[10]; 
out[19] += Buffer; 
 
Buffer = (-0.5652653816)*in[10]; 
out[20] += Buffer; 
 
Buffer = (-5.1517105306)*in[10]; 
out[23] += Buffer; 
 
Buffer = (4.1053953962)*in[10]; 
out[25] += Buffer; 
 
Buffer = (-1.4719104438)*in[10]; 
out[26] += Buffer; 
 
Buffer = (4.1053953962)*in[10]; 
out[27] += Buffer; 
 
Buffer = (-1.4719104438)*in[10]; 
out[29] += Buffer; 
 
Buffer = (-0.4110279866)*in[10]; 
out[30] += Buffer; 
 
Buffer = (-0.4110279866)*in[10]; 
out[31] += Buffer; 
 
Buffer = (-0.0559058514)*in[10]; 
out[32] += Buffer; 
 
Buffer = (-4.1053953962)*in[10]; 
out[33] += Buffer; 
 
Buffer = (-4.1053953962)*in[10]; 
out[34] += Buffer; 
 
Buffer = (-5.5773058400)*in[10]; 
out[35] += Buffer; 
 
Buffer = (0.4110279866)*in[10]; 
out[37] += Buffer; 
 
Buffer = (0.3551221352)*in[10]; 
out[38] += Buffer; 
 
Buffer = (0.4110279866)*in[10]; 
out[39] += Buffer; 
 
Buffer = (0.3551221352)*in[10]; 
out[41] += Buffer; 
 
Buffer = (0.9908122592)*in[11]; 
out[0] += Buffer; 
 
Buffer = (-0.8164535405)*in[11]; 
out[1] += Buffer; 
 
Buffer = (0.9908122592)*in[11]; 
out[2] += Buffer; 
 
Buffer = (-0.9908122592)*in[11]; 
out[3] += Buffer; 
 
Buffer = (-1.8072657997)*in[11]; 
out[4] += Buffer; 
 
Buffer = (0.3520629014)*in[11]; 
out[7] += Buffer; 
 
Buffer = (-1.8072657997)*in[11]; 
out[10] += Buffer; 
 
Buffer = (-0.9908122592)*in[11]; 
out[11] += Buffer; 
 
Buffer = (0.4978729532)*in[11]; 
out[12] += Buffer; 
 
Buffer = (-0.0673924284)*in[11]; 
out[13] += Buffer; 
 
Buffer = (0.4978729532)*in[11]; 
out[14] += Buffer; 
 
Buffer = (-0.4978729532)*in[11]; 
out[15] += Buffer; 
 
Buffer = (-0.5652653816)*in[11]; 
out[16] += Buffer; 
 
Buffer = (-5.1517105306)*in[11]; 
out[19] += Buffer; 
 
Buffer = (-0.5652653816)*in[11]; 
out[22] += Buffer; 
 
Buffer = (-0.4978729532)*in[11]; 
out[23] += Buffer; 
 
Buffer = (-1.4719104438)*in[11]; 
out[25] += Buffer; 
 
Buffer = (4.1053953962)*in[11]; 
out[26] += Buffer; 
 
Buffer = (-0.4110279866)*in[11]; 
out[27] += Buffer; 
 
Buffer = (-0.0559058514)*in[11]; 
out[28] += Buffer; 
 
Buffer = (-0.4110279866)*in[11]; 
out[29] += Buffer; 
 
Buffer = (4.1053953962)*in[11]; 
out[30] += Buffer; 
 
Buffer = (-1.4719104438)*in[11]; 
out[31] += Buffer; 
 
Buffer = (0.4110279866)*in[11]; 
out[33] += Buffer; 
 
Buffer = (0.3551221352)*in[11]; 
out[34] += Buffer; 
 
Buffer = (0.3551221352)*in[11]; 
out[37] += Buffer; 
 
Buffer = (0.4110279866)*in[11]; 
out[38] += Buffer; 
 
Buffer = (-4.1053953962)*in[11]; 
out[39] += Buffer; 
 
Buffer = (-5.5773058400)*in[11]; 
out[40] += Buffer; 
 
Buffer = (-4.1053953962)*in[11]; 
out[41] += Buffer; 
 
Buffer = (-0.8164535405)*in[12]; 
out[0] += Buffer; 
 
Buffer = (0.9908122592)*in[12]; 
out[1] += Buffer; 
 
Buffer = (0.9908122592)*in[12]; 
out[2] += Buffer; 
 
Buffer = (0.3520629014)*in[12]; 
out[3] += Buffer; 
 
Buffer = (-1.8072657997)*in[12]; 
out[6] += Buffer; 
 
Buffer = (-0.9908122592)*in[12]; 
out[7] += Buffer; 
 
Buffer = (-1.8072657997)*in[12]; 
out[9] += Buffer; 
 
Buffer = (-0.9908122592)*in[12]; 
out[11] += Buffer; 
 
Buffer = (-0.0673924284)*in[12]; 
out[12] += Buffer; 
 
Buffer = (0.4978729532)*in[12]; 
out[13] += Buffer; 
 
Buffer = (0.4978729532)*in[12]; 
out[14] += Buffer; 
 
Buffer = (-5.1517105306)*in[12]; 
out[15] += Buffer; 
 
Buffer = (-0.5652653816)*in[12]; 
out[18] += Buffer; 
 
Buffer = (-0.4978729532)*in[12]; 
out[19] += Buffer; 
 
Buffer = (-0.5652653816)*in[12]; 
out[21] += Buffer; 
 
Buffer = (-0.4978729532)*in[12]; 
out[23] += Buffer; 
 
Buffer = (-0.0559058514)*in[12]; 
out[24] += Buffer; 
 
Buffer = (-0.4110279866)*in[12]; 
out[25] += Buffer; 
 
Buffer = (-0.4110279866)*in[12]; 
out[26] += Buffer; 
 
Buffer = (-1.4719104438)*in[12]; 
out[27] += Buffer; 
 
Buffer = (4.1053953962)*in[12]; 
out[29] += Buffer; 
 
Buffer = (-1.4719104438)*in[12]; 
out[30] += Buffer; 
 
Buffer = (4.1053953962)*in[12]; 
out[31] += Buffer; 
 
Buffer = (0.3551221352)*in[12]; 
out[33] += Buffer; 
 
Buffer = (0.4110279866)*in[12]; 
out[34] += Buffer; 
 
Buffer = (-5.5773058400)*in[12]; 
out[36] += Buffer; 
 
Buffer = (-4.1053953962)*in[12]; 
out[37] += Buffer; 
 
Buffer = (-4.1053953962)*in[12]; 
out[38] += Buffer; 
 
Buffer = (0.3551221352)*in[12]; 
out[39] += Buffer; 
 
Buffer = (0.4110279866)*in[12]; 
out[41] += Buffer; 
 
Buffer = (-0.3520629014)*in[13]; 
out[0] += Buffer; 
 
Buffer = (-0.3520629014)*in[13]; 
out[1] += Buffer; 
 
Buffer = (-0.3520629014)*in[13]; 
out[2] += Buffer; 
 
Buffer = (0.8164535405)*in[13]; 
out[3] += Buffer; 
 
Buffer = (1.8072657997)*in[13]; 
out[4] += Buffer; 
 
Buffer = (1.8072657997)*in[13]; 
out[5] += Buffer; 
 
Buffer = (1.8072657997)*in[13]; 
out[6] += Buffer; 
 
Buffer = (0.8164535405)*in[13]; 
out[7] += Buffer; 
 
Buffer = (1.8072657997)*in[13]; 
out[8] += Buffer; 
 
Buffer = (1.8072657997)*in[13]; 
out[9] += Buffer; 
 
Buffer = (1.8072657997)*in[13]; 
out[10] += Buffer; 
 
Buffer = (0.8164535405)*in[13]; 
out[11] += Buffer; 
 
Buffer = (5.1517105306)*in[13]; 
out[12] += Buffer; 
 
Buffer = (5.1517105306)*in[13]; 
out[13] += Buffer; 
 
Buffer = (5.1517105306)*in[13]; 
out[14] += Buffer; 
 
Buffer = (0.0673924284)*in[13]; 
out[15] += Buffer; 
 
Buffer = (0.5652653816)*in[13]; 
out[16] += Buffer; 
 
Buffer = (0.5652653816)*in[13]; 
out[17] += Buffer; 
 
Buffer = (0.5652653816)*in[13]; 
out[18] += Buffer; 
 
Buffer = (0.0673924284)*in[13]; 
out[19] += Buffer; 
 
Buffer = (0.5652653816)*in[13]; 
out[20] += Buffer; 
 
Buffer = (0.5652653816)*in[13]; 
out[21] += Buffer; 
 
Buffer = (0.5652653816)*in[13]; 
out[22] += Buffer; 
 
Buffer = (0.0673924284)*in[13]; 
out[23] += Buffer; 
 
Buffer = (0.0559058514)*in[13]; 
out[24] += Buffer; 
 
Buffer = (-0.3551221352)*in[13]; 
out[25] += Buffer; 
 
Buffer = (-0.3551221352)*in[13]; 
out[26] += Buffer; 
 
Buffer = (-0.3551221352)*in[13]; 
out[27] += Buffer; 
 
Buffer = (0.0559058514)*in[13]; 
out[28] += Buffer; 
 
Buffer = (-0.3551221352)*in[13]; 
out[29] += Buffer; 
 
Buffer = (-0.3551221352)*in[13]; 
out[30] += Buffer; 
 
Buffer = (-0.3551221352)*in[13]; 
out[31] += Buffer; 
 
Buffer = (0.0559058514)*in[13]; 
out[32] += Buffer; 
 
Buffer = (1.4719104438)*in[13]; 
out[33] += Buffer; 
 
Buffer = (1.4719104438)*in[13]; 
out[34] += Buffer; 
 
Buffer = (5.5773058400)*in[13]; 
out[35] += Buffer; 
 
Buffer = (5.5773058400)*in[13]; 
out[36] += Buffer; 
 
Buffer = (1.4719104438)*in[13]; 
out[37] += Buffer; 
 
Buffer = (1.4719104438)*in[13]; 
out[38] += Buffer; 
 
Buffer = (1.4719104438)*in[13]; 
out[39] += Buffer; 
 
Buffer = (5.5773058400)*in[13]; 
out[40] += Buffer; 
 
Buffer = (1.4719104438)*in[13]; 
out[41] += Buffer; 
 
Buffer = (1.3849083756)*in[14]; 
out[0] += Buffer; 
 
Buffer = (1.3849083756)*in[14]; 
out[1] += Buffer; 
 
Buffer = (1.3849083756)*in[14]; 
out[2] += Buffer; 
 
Buffer = (-1.3849083756)*in[14]; 
out[3] += Buffer; 
 
Buffer = (-1.3849083756)*in[14]; 
out[7] += Buffer; 
 
Buffer = (-1.3849083756)*in[14]; 
out[11] += Buffer; 
 
Buffer = (-6.0258584915)*in[14]; 
out[12] += Buffer; 
 
Buffer = (-6.0258584915)*in[14]; 
out[13] += Buffer; 
 
Buffer = (-6.0258584915)*in[14]; 
out[14] += Buffer; 
 
Buffer = (6.0258584915)*in[14]; 
out[15] += Buffer; 
 
Buffer = (6.0258584915)*in[14]; 
out[19] += Buffer; 
 
Buffer = (6.0258584915)*in[14]; 
out[23] += Buffer; 
 
Buffer = (2.1653719455)*in[14]; 
out[25] += Buffer; 
 
Buffer = (2.1653719455)*in[14]; 
out[26] += Buffer; 
 
Buffer = (2.1653719455)*in[14]; 
out[27] += Buffer; 
 
Buffer = (2.1653719455)*in[14]; 
out[29] += Buffer; 
 
Buffer = (2.1653719455)*in[14]; 
out[30] += Buffer; 
 
Buffer = (2.1653719455)*in[14]; 
out[31] += Buffer; 
 
Buffer = (-2.1653719455)*in[14]; 
out[33] += Buffer; 
 
Buffer = (-2.1653719455)*in[14]; 
out[34] += Buffer; 
 
Buffer = (-2.1653719455)*in[14]; 
out[37] += Buffer; 
 
Buffer = (-2.1653719455)*in[14]; 
out[38] += Buffer; 
 
Buffer = (-2.1653719455)*in[14]; 
out[39] += Buffer; 
 
Buffer = (-2.1653719455)*in[14]; 
out[41] += Buffer; 
 


 

  
            }

            inline void TransposeGradientInterpolation_Opt(const Real * in, Real *out)
            {
                
            }

            inline void GradientInterpolation_Opt(const Real * in, Real *out)
            {
                Real Buffer;

                for (Index i=0;i<42;++i) 
                {
                    out[i]=0.0;
                }


Buffer = (-1.5846000073282709)*in[0]; 
out[0] += Buffer; 
out[1] += Buffer; 
out[2] += Buffer; 

Buffer = (0.2748830499140588)*in[0]; 
out[3] += Buffer; 
out[7] += Buffer; 
out[11] += Buffer; 

Buffer = (0.4282516876757742)*in[0]; 
out[4] += Buffer; 
out[5] += Buffer; 
out[6] += Buffer; 
out[8] += Buffer; 
out[9] += Buffer; 
out[10] += Buffer; 

Buffer = (0.0805491061721774)*in[0]; 
out[12] += Buffer; 
out[13] += Buffer; 
out[14] += Buffer; 

Buffer = (-0.1733620808428973)*in[0]; 
out[15] += Buffer; 
out[19] += Buffer; 
out[23] += Buffer; 

Buffer = (-0.3063509416774171)*in[0]; 
out[16] += Buffer; 
out[17] += Buffer; 
out[18] += Buffer; 
out[20] += Buffer; 
out[21] += Buffer; 
out[22] += Buffer; 

Buffer = (-0.8241969447648654)*in[0]; 
out[24] += Buffer; 
out[28] += Buffer; 
out[32] += Buffer; 

Buffer = (-0.3003740274442528)*in[0]; 
out[25] += Buffer; 
out[26] += Buffer; 
out[27] += Buffer; 
out[29] += Buffer; 
out[30] += Buffer; 
out[31] += Buffer; 

Buffer = (0.1762864376022318)*in[0]; 
out[33] += Buffer; 
out[34] += Buffer; 
out[37] += Buffer; 
out[38] += Buffer; 
out[39] += Buffer; 
out[41] += Buffer; 

Buffer = (0.1982845346068864)*in[0]; 
out[35] += Buffer; 
out[36] += Buffer; 
out[40] += Buffer; 

Buffer = (-0.2748830499140588)*in[1]; 
out[0] += Buffer; 

Buffer = (0.1533686377617155)*in[1]; 
out[1] += Buffer; 
out[2] += Buffer; 
out[7] -= Buffer; 
out[11] -= Buffer; 

Buffer = (1.5846000073282709)*in[1]; 
out[3] += Buffer; 

Buffer = (-0.4282516876757742)*in[1]; 
out[6] += Buffer; 
out[9] += Buffer; 

Buffer = (0.1733620808428973)*in[1]; 
out[12] += Buffer; 

Buffer = (-0.1329888608345198)*in[1]; 
out[13] += Buffer; 
out[14] += Buffer; 
out[19] -= Buffer; 
out[23] -= Buffer; 

Buffer = (-0.0805491061721774)*in[1]; 
out[15] += Buffer; 

Buffer = (0.3063509416774171)*in[1]; 
out[18] += Buffer; 
out[21] += Buffer; 

Buffer = (0.8241969447648654)*in[1]; 
out[24] += Buffer; 

Buffer = (0.5238229173206126)*in[1]; 
out[25] += Buffer; 
out[26] += Buffer; 
out[34] -= Buffer; 
out[41] -= Buffer; 

Buffer = (-0.1762864376022318)*in[1]; 
out[27] += Buffer; 
out[30] += Buffer; 

Buffer = (0.0219980970046546)*in[1]; 
out[29] += Buffer; 
out[31] += Buffer; 
out[37] -= Buffer; 
out[38] -= Buffer; 

Buffer = (0.3003740274442528)*in[1]; 
out[33] += Buffer; 
out[39] += Buffer; 

Buffer = (-0.1982845346068864)*in[1]; 
out[36] += Buffer; 

Buffer = (0.1533686377617155)*in[2]; 
out[0] += Buffer; 
out[2] += Buffer; 
out[3] -= Buffer; 
out[11] -= Buffer; 

Buffer = (-0.2748830499140588)*in[2]; 
out[1] += Buffer; 

Buffer = (-0.4282516876757742)*in[2]; 
out[4] += Buffer; 
out[10] += Buffer; 

Buffer = (1.5846000073282709)*in[2]; 
out[7] += Buffer; 

Buffer = (-0.1329888608345198)*in[2]; 
out[12] += Buffer; 
out[14] += Buffer; 
out[15] -= Buffer; 
out[23] -= Buffer; 

Buffer = (0.1733620808428973)*in[2]; 
out[13] += Buffer; 

Buffer = (0.3063509416774171)*in[2]; 
out[16] += Buffer; 
out[22] += Buffer; 

Buffer = (-0.0805491061721774)*in[2]; 
out[19] += Buffer; 

Buffer = (-0.1762864376022318)*in[2]; 
out[25] += Buffer; 
out[31] += Buffer; 

Buffer = (0.0219980970046546)*in[2]; 
out[26] += Buffer; 
out[30] += Buffer; 
out[39] -= Buffer; 
out[41] -= Buffer; 

Buffer = (0.5238229173206126)*in[2]; 
out[27] += Buffer; 
out[29] += Buffer; 
out[33] -= Buffer; 
out[38] -= Buffer; 

Buffer = (0.8241969447648654)*in[2]; 
out[28] += Buffer; 

Buffer = (0.3003740274442528)*in[2]; 
out[34] += Buffer; 
out[37] += Buffer; 

Buffer = (-0.1982845346068864)*in[2]; 
out[40] += Buffer; 

Buffer = (0.1533686377617155)*in[3]; 
out[0] += Buffer; 
out[1] += Buffer; 
out[3] -= Buffer; 
out[7] -= Buffer; 

Buffer = (-0.2748830499140588)*in[3]; 
out[2] += Buffer; 

Buffer = (-0.4282516876757742)*in[3]; 
out[5] += Buffer; 
out[8] += Buffer; 

Buffer = (1.5846000073282709)*in[3]; 
out[11] += Buffer; 

Buffer = (-0.1329888608345198)*in[3]; 
out[12] += Buffer; 
out[13] += Buffer; 
out[15] -= Buffer; 
out[19] -= Buffer; 

Buffer = (0.1733620808428973)*in[3]; 
out[14] += Buffer; 

Buffer = (0.3063509416774171)*in[3]; 
out[17] += Buffer; 
out[20] += Buffer; 

Buffer = (-0.0805491061721774)*in[3]; 
out[23] += Buffer; 

Buffer = (0.0219980970046546)*in[3]; 
out[25] += Buffer; 
out[27] += Buffer; 
out[33] -= Buffer; 
out[34] -= Buffer; 

Buffer = (-0.1762864376022318)*in[3]; 
out[26] += Buffer; 
out[29] += Buffer; 

Buffer = (0.5238229173206126)*in[3]; 
out[30] += Buffer; 
out[31] += Buffer; 
out[37] -= Buffer; 
out[39] -= Buffer; 

Buffer = (0.8241969447648654)*in[3]; 
out[32] += Buffer; 

Buffer = (-0.1982845346068864)*in[3]; 
out[35] += Buffer; 

Buffer = (0.3003740274442528)*in[3]; 
out[38] += Buffer; 
out[41] += Buffer; 

Buffer = (1.2892868929320165)*in[4]; 
out[0] += Buffer; 
out[3] -= Buffer; 

Buffer = (-0.7946608590146454)*in[4]; 
out[1] += Buffer; 
out[2] += Buffer; 

Buffer = (-2.0839477519466620)*in[4]; 
out[4] += Buffer; 
out[5] += Buffer; 

Buffer = (0.2591746961328832)*in[4]; 
out[7] += Buffer; 
out[11] += Buffer; 

Buffer = (0.4322882430790784)*in[4]; 
out[8] += Buffer; 
out[10] += Buffer; 

Buffer = (0.0897361784633483)*in[4]; 
out[12] += Buffer; 
out[15] -= Buffer; 

Buffer = (0.0715962681198142)*in[4]; 
out[13] += Buffer; 
out[14] += Buffer; 

Buffer = (-0.0181399103435341)*in[4]; 
out[16] += Buffer; 
out[17] += Buffer; 

Buffer = (-0.2390823071177199)*in[4]; 
out[19] += Buffer; 
out[23] += Buffer; 

Buffer = (-0.9923146185563440)*in[4]; 
out[20] += Buffer; 
out[22] += Buffer; 

Buffer = (-3.5297659265002421)*in[4]; 
out[25] += Buffer; 
out[26] += Buffer; 

Buffer = (-0.5473136903062660)*in[4]; 
out[27] += Buffer; 
out[30] += Buffer; 
out[33] -= Buffer; 
out[39] -= Buffer; 

Buffer = (-0.1571677714327426)*in[4]; 
out[28] += Buffer; 
out[32] += Buffer; 

Buffer = (0.1135037217583928)*in[4]; 
out[29] += Buffer; 
out[31] += Buffer; 

Buffer = (0.3901459188735233)*in[4]; 
out[34] += Buffer; 
out[41] += Buffer; 

Buffer = (0.6608174120646588)*in[4]; 
out[35] += Buffer; 
out[40] += Buffer; 

Buffer = (0.0188170046044938)*in[4]; 
out[37] += Buffer; 
out[38] += Buffer; 

Buffer = (-0.2591746961328832)*in[5]; 
out[0] += Buffer; 
out[1] += Buffer; 

Buffer = (0.1731135469461951)*in[5]; 
out[2] += Buffer; 
out[11] -= Buffer; 

Buffer = (0.7946608590146454)*in[5]; 
out[3] += Buffer; 
out[7] += Buffer; 

Buffer = (2.0839477519466620)*in[5]; 
out[4] += Buffer; 
out[6] += Buffer; 

Buffer = (-0.4322882430790784)*in[5]; 
out[9] += Buffer; 
out[10] += Buffer; 

Buffer = (0.2390823071177199)*in[5]; 
out[12] += Buffer; 
out[13] += Buffer; 

Buffer = (-0.7532323114386240)*in[5]; 
out[14] += Buffer; 
out[23] -= Buffer; 

Buffer = (-0.0715962681198142)*in[5]; 
out[15] += Buffer; 
out[19] += Buffer; 

Buffer = (0.0181399103435341)*in[5]; 
out[16] += Buffer; 
out[18] += Buffer; 

Buffer = (0.9923146185563440)*in[5]; 
out[21] += Buffer; 
out[22] += Buffer; 

Buffer = (0.1571677714327426)*in[5]; 
out[24] += Buffer; 
out[28] += Buffer; 

Buffer = (-0.3901459188735233)*in[5]; 
out[25] += Buffer; 
out[27] += Buffer; 

Buffer = (0.2706714931911355)*in[5]; 
out[26] += Buffer; 
out[29] += Buffer; 
out[38] -= Buffer; 
out[41] -= Buffer; 

Buffer = (-0.0188170046044938)*in[5]; 
out[30] += Buffer; 
out[31] += Buffer; 

Buffer = (3.5297659265002421)*in[5]; 
out[33] += Buffer; 
out[34] += Buffer; 

Buffer = (-0.6608174120646588)*in[5]; 
out[36] += Buffer; 
out[40] += Buffer; 

Buffer = (-0.1135037217583928)*in[5]; 
out[37] += Buffer; 
out[39] += Buffer; 

Buffer = (-0.7946608590146454)*in[6]; 
out[0] += Buffer; 
out[2] += Buffer; 

Buffer = (1.2892868929320165)*in[6]; 
out[1] += Buffer; 
out[7] -= Buffer; 

Buffer = (0.2591746961328832)*in[6]; 
out[3] += Buffer; 
out[11] += Buffer; 

Buffer = (0.4322882430790784)*in[6]; 
out[5] += Buffer; 
out[9] += Buffer; 

Buffer = (-2.0839477519466620)*in[6]; 
out[6] += Buffer; 
out[8] += Buffer; 

Buffer = (0.0715962681198142)*in[6]; 
out[12] += Buffer; 
out[14] += Buffer; 

Buffer = (0.0897361784633483)*in[6]; 
out[13] += Buffer; 
out[19] -= Buffer; 

Buffer = (-0.2390823071177199)*in[6]; 
out[15] += Buffer; 
out[23] += Buffer; 

Buffer = (-0.9923146185563440)*in[6]; 
out[17] += Buffer; 
out[21] += Buffer; 

Buffer = (-0.0181399103435341)*in[6]; 
out[18] += Buffer; 
out[20] += Buffer; 

Buffer = (-0.1571677714327426)*in[6]; 
out[24] += Buffer; 
out[32] += Buffer; 

Buffer = (-0.5473136903062660)*in[6]; 
out[25] += Buffer; 
out[31] += Buffer; 
out[34] -= Buffer; 
out[37] -= Buffer; 

Buffer = (0.1135037217583928)*in[6]; 
out[26] += Buffer; 
out[30] += Buffer; 

Buffer = (-3.5297659265002421)*in[6]; 
out[27] += Buffer; 
out[29] += Buffer; 

Buffer = (0.3901459188735233)*in[6]; 
out[33] += Buffer; 
out[38] += Buffer; 

Buffer = (0.6608174120646588)*in[6]; 
out[35] += Buffer; 
out[36] += Buffer; 

Buffer = (0.0188170046044938)*in[6]; 
out[39] += Buffer; 
out[41] += Buffer; 

Buffer = (-0.7946608590146454)*in[7]; 
out[0] += Buffer; 
out[1] += Buffer; 

Buffer = (1.2892868929320165)*in[7]; 
out[2] += Buffer; 
out[11] -= Buffer; 

Buffer = (0.2591746961328832)*in[7]; 
out[3] += Buffer; 
out[7] += Buffer; 

Buffer = (0.4322882430790784)*in[7]; 
out[4] += Buffer; 
out[6] += Buffer; 

Buffer = (-2.0839477519466620)*in[7]; 
out[9] += Buffer; 
out[10] += Buffer; 

Buffer = (0.0715962681198142)*in[7]; 
out[12] += Buffer; 
out[13] += Buffer; 

Buffer = (0.0897361784633483)*in[7]; 
out[14] += Buffer; 
out[23] -= Buffer; 

Buffer = (-0.2390823071177199)*in[7]; 
out[15] += Buffer; 
out[19] += Buffer; 

Buffer = (-0.9923146185563440)*in[7]; 
out[16] += Buffer; 
out[18] += Buffer; 

Buffer = (-0.0181399103435341)*in[7]; 
out[21] += Buffer; 
out[22] += Buffer; 

Buffer = (-0.1571677714327426)*in[7]; 
out[24] += Buffer; 
out[28] += Buffer; 

Buffer = (0.1135037217583928)*in[7]; 
out[25] += Buffer; 
out[27] += Buffer; 

Buffer = (-0.5473136903062660)*in[7]; 
out[26] += Buffer; 
out[29] += Buffer; 
out[38] -= Buffer; 
out[41] -= Buffer; 

Buffer = (-3.5297659265002421)*in[7]; 
out[30] += Buffer; 
out[31] += Buffer; 

Buffer = (0.0188170046044938)*in[7]; 
out[33] += Buffer; 
out[34] += Buffer; 

Buffer = (0.6608174120646588)*in[7]; 
out[36] += Buffer; 
out[40] += Buffer; 

Buffer = (0.3901459188735233)*in[7]; 
out[37] += Buffer; 
out[39] += Buffer; 

Buffer = (-0.2591746961328832)*in[8]; 
out[0] += Buffer; 
out[2] += Buffer; 

Buffer = (0.1731135469461951)*in[8]; 
out[1] += Buffer; 
out[7] -= Buffer; 

Buffer = (0.7946608590146454)*in[8]; 
out[3] += Buffer; 
out[11] += Buffer; 

Buffer = (2.0839477519466620)*in[8]; 
out[5] += Buffer; 
out[9] += Buffer; 

Buffer = (-0.4322882430790784)*in[8]; 
out[6] += Buffer; 
out[8] += Buffer; 

Buffer = (0.2390823071177199)*in[8]; 
out[12] += Buffer; 
out[14] += Buffer; 

Buffer = (-0.7532323114386240)*in[8]; 
out[13] += Buffer; 
out[19] -= Buffer; 

Buffer = (-0.0715962681198142)*in[8]; 
out[15] += Buffer; 
out[23] += Buffer; 

Buffer = (0.0181399103435341)*in[8]; 
out[17] += Buffer; 
out[21] += Buffer; 

Buffer = (0.9923146185563440)*in[8]; 
out[18] += Buffer; 
out[20] += Buffer; 

Buffer = (0.1571677714327426)*in[8]; 
out[24] += Buffer; 
out[32] += Buffer; 

Buffer = (0.2706714931911355)*in[8]; 
out[25] += Buffer; 
out[31] += Buffer; 
out[34] -= Buffer; 
out[37] -= Buffer; 

Buffer = (-0.3901459188735233)*in[8]; 
out[26] += Buffer; 
out[30] += Buffer; 

Buffer = (-0.0188170046044938)*in[8]; 
out[27] += Buffer; 
out[29] += Buffer; 

Buffer = (-0.1135037217583928)*in[8]; 
out[33] += Buffer; 
out[38] += Buffer; 

Buffer = (-0.6608174120646588)*in[8]; 
out[35] += Buffer; 
out[36] += Buffer; 

Buffer = (3.5297659265002421)*in[8]; 
out[39] += Buffer; 
out[41] += Buffer; 

Buffer = (0.1731135469461951)*in[9]; 
out[0] += Buffer; 
out[3] -= Buffer; 

Buffer = (-0.2591746961328832)*in[9]; 
out[1] += Buffer; 
out[2] += Buffer; 

Buffer = (-0.4322882430790784)*in[9]; 
out[4] += Buffer; 
out[5] += Buffer; 

Buffer = (0.7946608590146454)*in[9]; 
out[7] += Buffer; 
out[11] += Buffer; 

Buffer = (2.0839477519466620)*in[9]; 
out[8] += Buffer; 
out[10] += Buffer; 

Buffer = (-0.7532323114386240)*in[9]; 
out[12] += Buffer; 
out[15] -= Buffer; 

Buffer = (0.2390823071177199)*in[9]; 
out[13] += Buffer; 
out[14] += Buffer; 

Buffer = (0.9923146185563440)*in[9]; 
out[16] += Buffer; 
out[17] += Buffer; 

Buffer = (-0.0715962681198142)*in[9]; 
out[19] += Buffer; 
out[23] += Buffer; 

Buffer = (0.0181399103435341)*in[9]; 
out[20] += Buffer; 
out[22] += Buffer; 

Buffer = (-0.0188170046044938)*in[9]; 
out[25] += Buffer; 
out[26] += Buffer; 

Buffer = (0.2706714931911355)*in[9]; 
out[27] += Buffer; 
out[30] += Buffer; 
out[33] -= Buffer; 
out[39] -= Buffer; 

Buffer = (0.1571677714327426)*in[9]; 
out[28] += Buffer; 
out[32] += Buffer; 

Buffer = (-0.3901459188735233)*in[9]; 
out[29] += Buffer; 
out[31] += Buffer; 

Buffer = (-0.1135037217583928)*in[9]; 
out[34] += Buffer; 
out[41] += Buffer; 

Buffer = (-0.6608174120646588)*in[9]; 
out[35] += Buffer; 
out[40] += Buffer; 

Buffer = (3.5297659265002421)*in[9]; 
out[37] += Buffer; 
out[38] += Buffer; 

Buffer = (0.9908122592265003)*in[10]; 
out[0] += Buffer; 
out[1] += Buffer; 
out[3] -= Buffer; 
out[7] -= Buffer; 

Buffer = (-0.8164535404994472)*in[10]; 
out[2] += Buffer; 

Buffer = (-1.8072657997259474)*in[10]; 
out[5] += Buffer; 
out[8] += Buffer; 

Buffer = (0.3520629013873700)*in[10]; 
out[11] += Buffer; 

Buffer = (0.4978729532262263)*in[10]; 
out[12] += Buffer; 
out[13] += Buffer; 
out[15] -= Buffer; 
out[19] -= Buffer; 

Buffer = (-0.0673924283917054)*in[10]; 
out[14] += Buffer; 

Buffer = (-0.5652653816179317)*in[10]; 
out[17] += Buffer; 
out[20] += Buffer; 

Buffer = (-5.1517105306024176)*in[10]; 
out[23] += Buffer; 

Buffer = (4.1053953962054592)*in[10]; 
out[25] += Buffer; 
out[27] += Buffer; 
out[33] -= Buffer; 
out[34] -= Buffer; 

Buffer = (-1.4719104438091766)*in[10]; 
out[26] += Buffer; 
out[29] += Buffer; 

Buffer = (-0.4110279866381638)*in[10]; 
out[30] += Buffer; 
out[31] += Buffer; 
out[37] -= Buffer; 
out[39] -= Buffer; 

Buffer = (-0.0559058514071760)*in[10]; 
out[32] += Buffer; 

Buffer = (-5.5773058400146356)*in[10]; 
out[35] += Buffer; 

Buffer = (0.3551221352309878)*in[10]; 
out[38] += Buffer; 
out[41] += Buffer; 

Buffer = (0.9908122592265003)*in[11]; 
out[0] += Buffer; 
out[2] += Buffer; 
out[3] -= Buffer; 
out[11] -= Buffer; 

Buffer = (-0.8164535404994472)*in[11]; 
out[1] += Buffer; 

Buffer = (-1.8072657997259474)*in[11]; 
out[4] += Buffer; 
out[10] += Buffer; 

Buffer = (0.3520629013873700)*in[11]; 
out[7] += Buffer; 

Buffer = (0.4978729532262263)*in[11]; 
out[12] += Buffer; 
out[14] += Buffer; 
out[15] -= Buffer; 
out[23] -= Buffer; 

Buffer = (-0.0673924283917054)*in[11]; 
out[13] += Buffer; 

Buffer = (-0.5652653816179317)*in[11]; 
out[16] += Buffer; 
out[22] += Buffer; 

Buffer = (-5.1517105306024176)*in[11]; 
out[19] += Buffer; 

Buffer = (-1.4719104438091766)*in[11]; 
out[25] += Buffer; 
out[31] += Buffer; 

Buffer = (4.1053953962054592)*in[11]; 
out[26] += Buffer; 
out[30] += Buffer; 
out[39] -= Buffer; 
out[41] -= Buffer; 

Buffer = (-0.4110279866381638)*in[11]; 
out[27] += Buffer; 
out[29] += Buffer; 
out[33] -= Buffer; 
out[38] -= Buffer; 

Buffer = (-0.0559058514071760)*in[11]; 
out[28] += Buffer; 

Buffer = (0.3551221352309878)*in[11]; 
out[34] += Buffer; 
out[37] += Buffer; 

Buffer = (-5.5773058400146356)*in[11]; 
out[40] += Buffer; 

Buffer = (-0.8164535404994472)*in[12]; 
out[0] += Buffer; 

Buffer = (0.9908122592265003)*in[12]; 
out[1] += Buffer; 
out[2] += Buffer; 
out[7] -= Buffer; 
out[11] -= Buffer; 

Buffer = (0.3520629013873700)*in[12]; 
out[3] += Buffer; 

Buffer = (-1.8072657997259474)*in[12]; 
out[6] += Buffer; 
out[9] += Buffer; 

Buffer = (-0.0673924283917054)*in[12]; 
out[12] += Buffer; 

Buffer = (0.4978729532262263)*in[12]; 
out[13] += Buffer; 
out[14] += Buffer; 
out[19] -= Buffer; 
out[23] -= Buffer; 

Buffer = (-5.1517105306024176)*in[12]; 
out[15] += Buffer; 

Buffer = (-0.5652653816179317)*in[12]; 
out[18] += Buffer; 
out[21] += Buffer; 

Buffer = (-0.0559058514071760)*in[12]; 
out[24] += Buffer; 

Buffer = (-0.4110279866381638)*in[12]; 
out[25] += Buffer; 
out[26] += Buffer; 
out[34] -= Buffer; 
out[41] -= Buffer; 

Buffer = (-1.4719104438091766)*in[12]; 
out[27] += Buffer; 
out[30] += Buffer; 

Buffer = (4.1053953962054592)*in[12]; 
out[29] += Buffer; 
out[31] += Buffer; 
out[37] -= Buffer; 
out[38] -= Buffer; 

Buffer = (0.3551221352309878)*in[12]; 
out[33] += Buffer; 
out[39] += Buffer; 

Buffer = (-5.5773058400146356)*in[12]; 
out[36] += Buffer; 

Buffer = (-0.3520629013873700)*in[13]; 
out[0] += Buffer; 
out[1] += Buffer; 
out[2] += Buffer; 

Buffer = (0.8164535404994472)*in[13]; 
out[3] += Buffer; 
out[7] += Buffer; 
out[11] += Buffer; 

Buffer = (1.8072657997259474)*in[13]; 
out[4] += Buffer; 
out[5] += Buffer; 
out[6] += Buffer; 
out[8] += Buffer; 
out[9] += Buffer; 
out[10] += Buffer; 

Buffer = (5.1517105306024176)*in[13]; 
out[12] += Buffer; 
out[13] += Buffer; 
out[14] += Buffer; 

Buffer = (0.0673924283917054)*in[13]; 
out[15] += Buffer; 
out[19] += Buffer; 
out[23] += Buffer; 

Buffer = (0.5652653816179317)*in[13]; 
out[16] += Buffer; 
out[17] += Buffer; 
out[18] += Buffer; 
out[20] += Buffer; 
out[21] += Buffer; 
out[22] += Buffer; 

Buffer = (0.0559058514071760)*in[13]; 
out[24] += Buffer; 
out[28] += Buffer; 
out[32] += Buffer; 

Buffer = (-0.3551221352309878)*in[13]; 
out[25] += Buffer; 
out[26] += Buffer; 
out[27] += Buffer; 
out[29] += Buffer; 
out[30] += Buffer; 
out[31] += Buffer; 

Buffer = (1.4719104438091766)*in[13]; 
out[33] += Buffer; 
out[34] += Buffer; 
out[37] += Buffer; 
out[38] += Buffer; 
out[39] += Buffer; 
out[41] += Buffer; 

Buffer = (5.5773058400146356)*in[13]; 
out[35] += Buffer; 
out[36] += Buffer; 
out[40] += Buffer; 

Buffer = (1.3849083755695610)*in[14]; 
out[0] += Buffer; 
out[1] += Buffer; 
out[2] += Buffer; 
out[3] -= Buffer; 
out[7] -= Buffer; 
out[11] -= Buffer; 

Buffer = (-6.0258584915089921)*in[14]; 
out[12] += Buffer; 
out[13] += Buffer; 
out[14] += Buffer; 
out[15] -= Buffer; 
out[19] -= Buffer; 
out[23] -= Buffer; 

Buffer = (2.1653719455290839)*in[14]; 
out[25] += Buffer; 
out[26] += Buffer; 
out[27] += Buffer; 
out[29] += Buffer; 
out[30] += Buffer; 
out[31] += Buffer; 
out[33] -= Buffer; 
out[34] -= Buffer; 
out[37] -= Buffer; 
out[38] -= Buffer; 
out[39] -= Buffer; 
out[41] -= Buffer; 

            }

        }
    }
        
}
// -----------------------------------------------------------------------------------------//

 