
#include "../../paul.h"

double get_eint( double * , double * );
void get_Ederivs( double * , double * , double * );

struct NSE{
  double rho;  //density g/cm^3
  double T;
  double e;   //specific internal energy
  double etot;  //specific internal energy + specific binding energy
  double xn[NUM_I];
  double xn_exp[NUM_I];
  double yn[NUM_I];
  double log_data[NUM_I];
  double dlogX_dT9[NUM_I];

  double y_e;
  double y_e_init;
  double mu_p;
  double mu_n;
};

static int npts_1 = 72;

//Modes - TODO add to in.par
int solve_nse_e_mode = 1;
// 0 for (rho, T, Ye) or 1 for (rho, e, Ye) //
int input = 1; 

// this is T9
double temp_array[] = {
  0.01, 0.15, 0.2, 0.3, 0.4,
  0.5, 0.6, 0.7, 0.8, 0.9,
  1.0, 1.5, 2.0, 2.5, 3.0,
  3.5, 4.0, 4.5, 5.0, 6.0,
  7.0, 8.0, 9.0, 10.0, 12.0,
  14.0, 16.0, 18.0, 20.0, 22.0,
  24.0, 26.0, 28.0, 30.0, 35.0,
  40.0, 45.0, 50.0, 55.0, 60.0,
  65.0, 70.0, 75.0, 80.0, 85.0,
  90.0, 95.0, 100.0, 105.0, 110.0,
  115.0, 120.0, 125.0, 130.0, 135.0,
  140.0, 145.0, 150.0, 155.0, 160.0,
  165.0, 170.0, 175.0, 180.0, 190.0,
  200.0, 210.0, 220.0, 230.0, 240.0,
  250.0, 275.0,
};

// this is log10(partition function)
double O16_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.012837224705172217,
  0.03742649794062367, 0.07188200730612536, 0.1205739312058499, 0.1846914308175988, 0.26245108973042947,
  0.3463529744506387, 0.437750562820388, 0.534026106056135, 0.6344772701607315, 0.8981764834976765,
  1.1760912590556813, 1.4668676203541096, 1.7641761323903307, 2.0644579892269186, 2.367355921026019,
  2.667452952889954, 2.9656719712201065, 3.2624510897304293, 3.555094448578319, 3.845098040014257,
  4.133538908370218, 4.419955748489758, 4.704150516839799, 4.986771734266245, 5.267171728403014,
  5.547774705387822, 5.8267225201689925, 6.103803720955957, 6.380211241711606, 6.6551384348113825,
  6.929929560084588, 7.204119982655925, 7.477121254719663, 7.748962861256161, 8.021189299069938,
  8.292256071356476, 8.562292864456476, 8.832508912706237, 9.100370545117563, 9.640481436970422,
  10.178976947293169, 10.714329759745233, 11.250420002308894, 11.785329835010767, 12.320146286111054,
  12.856124444242301, 14.195899652409233,
};

// this is log10(partition function)
double Ne20_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 6.9486561213582446e-06, 0.00016586881316040883, 0.0011034421778731533, 0.003892457497077877,
  0.00954097493969645, 0.01859524021829981, 0.03107544483336982, 0.04661767038571622, 0.0846241727916796,
  0.12822183093465686, 0.174311933665943, 0.22124805254602342, 0.2683385291343481, 0.36172783601759284,
  0.456366033129043, 0.5514499979728752, 0.6483600109809317, 0.7466341989375788, 0.8481891169913987,
  0.9532763366673044, 1.0644579892269186, 1.1789769472931695, 1.3031960574204888, 1.6434526764861874,
  2.0170333392987803, 2.4099331233312946, 2.8068580295188172, 3.2013971243204513, 3.5899496013257077,
  3.9731278535996988, 4.352182518111363, 4.725911632295048, 5.096910013008056, 5.465382851448418,
  5.830588668685144, 6.193124598354461, 6.556302500767288, 6.916980047320382, 7.276461804173244,
  7.6344772701607315, 7.991669007379948, 8.34830486304816, 8.703291378118662, 9.056904851336473,
  9.411619705963231, 9.763427993562937, 10.117271295655764, 10.46686762035411, 10.818225893613956,
  11.170261715394957, 11.519827993775719, 11.869231719730976, 12.217483944213907, 12.916453948549925,
  13.613841821876068, 14.3096301674259, 15.004321373782643, 15.702430536445526, 16.399673721481037,
  17.096910013008056, 18.838849090737256,
};

// this is log10(partition function)
double Mg24_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  4.342942647204277e-07, 5.471765757979972e-05, 0.0007714899373308071, 0.0037633124724497633, 0.010764115210255056,
  0.022625058328435317, 0.039160607597355665, 0.05951911533271758, 0.08262238957783377, 0.13324118689139802,
  0.185518640557017, 0.2370005304649223, 0.2870228837145503, 0.3357157930198095, 0.43136376415898736,
  0.5263392773898441, 0.6253124509616739, 0.7307822756663892, 0.8463371121298052, 0.9749719942980689,
  1.1172712956557642, 1.2741578492636798, 1.4424797690644486, 1.6232492903979006, 2.103803720955957,
  2.598790506763115, 3.089905111439398, 3.577491799837225, 4.05307844348342, 4.52244423350632,
  4.984527313343793, 5.440909082065217, 5.894869656745253, 6.344392273685111, 6.791690649020118,
  7.235528446907549, 7.6785183790401135, 8.12057393120585, 8.558708570533165, 8.99563519459755,
  9.431363764158988, 9.866287339084195, 10.301029995663981, 10.732393759822969, 11.164352855784436,
  11.594392550375426, 12.02530586526477, 12.453318340047037, 12.881384656770573, 13.3096301674259,
  13.736396502276643, 14.161368002234974, 14.588831725594208, 15.012837224705173, 15.86569605991607,
  16.715167357848458, 17.56466606425209, 18.413299764081252, 19.26245108973043, 20.113943352306837,
  20.96284268120124, 23.089905111439396,
};

// this is log10(partition function)
double Si28_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 2.1714669808675565e-06, 7.12184552784347e-05, 0.0005624812393818785, 0.002223099674110693,
  0.0059171580771474625, 0.012282407118825528, 0.02157709561709228, 0.03370716078346824, 0.06502557053071237,
  0.1027522772573885, 0.14387160800291654, 0.1866035043986153, 0.2302807913268337, 0.3222192947339193,
  0.42324587393680785, 0.541579243946581, 0.6839471307515121, 0.8518696007297664, 1.041392685158225,
  1.250420002308894, 1.4727564493172123, 1.7024305364455252, 1.9375178920173466, 2.531478917042255,
  3.12057393120585, 3.7024305364455254, 4.271841606536499, 4.834420703681532, 5.389166084364533,
  5.937517892017347, 6.481442628502305, 7.021189299069938, 7.557507201905658, 8.089905111439398,
  8.622214022966295, 9.14921911265538, 9.675778341674086, 10.198657086954423, 10.721810615212547,
  11.2405492482826, 11.75966784468963, 12.276461804173245, 12.791690649020119, 13.305351369446624,
  13.818225893613956, 14.330413773349191, 14.840733234611807, 15.350248018334163, 15.85913829729453,
  16.367355921026018, 16.8750612633917, 17.38201704257487, 17.88874096068289, 18.90036712865647,
  19.911157608739977, 20.921166050637737, 21.9304395947667, 22.93951925261862, 23.948901760970212,
  24.958563883221967, 27.48572142648158,
};

// this is log10(partition function)
double S32_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 5.211502513843472e-06, 6.948155872801059e-05, 0.0003893875360542875,
  0.001336870159627728, 0.003378232401258555, 0.006963377556787149, 0.012456734172197398, 0.03011415790845077,
  0.057484285853877215, 0.0950053699501746, 0.14295136988131382, 0.20165707691270435, 0.3521825181113625,
  0.5502283530550941, 0.787460474518415, 1.0569048513364727, 1.3404441148401183, 1.631443769013172,
  1.92272545799326, 2.2121876044039577, 2.4955443375464483, 2.7737864449811935, 3.44870631990508,
  4.096910013008056, 4.726727209026572, 5.3404441148401185, 5.944975908412048, 6.541579243946581,
  7.133538908370218, 7.720985744153739, 8.305351369446624, 8.888740960682892, 9.469822015978163,
  10.049218022670182, 10.628388930050312, 11.20682587603185, 11.78175537465247, 12.356025857193123,
  12.9304395947667, 13.502427119984432, 14.075546961392531, 14.645422269349092, 15.214843848047698,
  15.783903579272735, 16.352182518111363, 16.920123326290724, 17.487138375477187, 18.053078443483418,
  18.621176281775035, 19.187520720836464, 19.753583058892907, 20.318063334962762, 21.450249108319362,
  22.580924975675618, 23.71264970162721, 24.84385542262316, 25.976808337338067, 27.110589710299248,
  28.24551266781415, 31.08635983067475,
};

// this is log10(partition function)
double Ar36_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 4.342942647204277e-07, 2.3451268844214655e-05, 0.00023141729162330258, 0.0010622869460975197,
  0.0031540913067783544, 0.007135153007315866, 0.013474284663478431, 0.02245187936733961, 0.048771089883939175,
  0.08643600351808534, 0.13560900039779808, 0.1965840257248699, 0.2696980636423851, 0.45331834004703764,
  0.6848453616444125, 0.9585638832219674, 1.2624510897304295, 1.5809249756756194, 1.9057958803678685,
  2.230448921378274, 2.550228353055094, 2.8662873390841948, 3.1760912590556813, 3.929418925714293,
  4.657055852857104, 5.365487984890899, 6.060697840353612, 6.746634198937579, 7.426511261364575,
  8.100370545117563, 8.773054693364262, 9.442479769064448, 10.11058971029925, 10.77451696572855,
  11.437750562820387, 12.100370545117563, 12.758911892397974, 13.41664050733828, 14.071882007306126,
  14.727541257028557, 15.38201704257487, 16.03342375548695, 16.684845361644413, 17.33445375115093,
  17.983626287124533, 18.63144376901317, 19.27875360095283, 19.92582757462474, 20.57170883180869,
  21.217483944213907, 21.863322860120455, 22.50785587169583, 23.152288344383056, 24.440909082065218,
  25.72916478969277, 27.01703333929878, 28.305351369446623, 29.595496221825574, 30.885926339801433,
  32.17897694729317, 35.41329976408125,
};

// this is log10(partition function)
double Ca40_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 3.4743419578801875e-06,
  2.6056887215373325e-05, 0.00012419046343446514, 0.0004254001802063995, 0.0011532564515138494, 0.005324252203746658,
  0.016451245325404363, 0.03938040551055626, 0.07909980819723089, 0.13972800117379408, 0.33041377334919086,
  0.6063813651106049, 0.9385197251764918, 1.2988530764097066, 1.6693168805661123, 2.037426497940624,
  2.403120521175818, 2.761927838420529, 3.113943352306837, 3.459392487759231, 4.301029995663981,
  5.117271295655764, 5.9148718175400505, 6.701567985055927, 7.478566495593843, 8.250420002308894,
  9.01703333929878, 9.781036938621131, 10.540329474790873, 11.296665190261532, 12.049218022670182,
  12.801403710017356, 13.549003262025789, 14.294466226161592, 15.037426497940624, 15.779596491257825,
  16.518513939877888, 17.255272505103306, 17.99211148778695, 18.72591163229505, 19.45939248775923,
  20.19033169817029, 20.920645001406786, 21.650307523131936, 22.378397900948137, 23.10720996964787,
  23.832508912706235, 24.558708570533167, 25.285557309007775, 26.008600171761916, 27.45939248775923,
  28.907948521612273, 30.356025857193124, 31.804820678721164, 33.25285303097989, 34.704150516839796,
  36.15533603746506, 39.78816837114117,
};

// this is log10(partition function)
double Ti44_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 4.342942647204277e-07, 1.7371744532199383e-06,
  7.382943437485089e-06, 0.0004987179011085028, 0.004043078170724821, 0.01413521502778782, 0.032426549056877405,
  0.05856115101668825, 0.09131586357749837, 0.1294359425571275, 0.17190802974603506, 0.2667731684215763,
  0.37035022176288673, 0.47788465213962983, 0.5860935485551829, 0.693748838923791, 0.9116901587538612,
  1.1522883443830565, 1.4409090820652177, 1.7788744720027396, 2.1522883443830563, 2.5415792439465807,
  2.9334872878487053, 3.322219294733919, 3.7041505168397992, 4.079181246047625, 4.996073654485276,
  5.885361220031512, 6.757396028793024, 7.619093330626742, 8.472756449317213, 9.32221929473392,
  10.167317334748176, 11.008600171761918, 11.85003325768977, 12.687528961214634, 13.52244423350632,
  14.354108439147401, 15.1846914308176, 16.012837224705173, 16.836956737059552, 17.65991620006985,
  18.481442628502304, 19.298853076409706, 20.117271295655765, 20.9329808219232, 21.746634198937578,
  22.559906625036113, 23.371067862271737, 24.181843587944773, 24.991226075692495, 25.799340549453582,
  26.60745502321467, 27.414973347970818, 28.220108088040057, 29.02530586526477, 30.63748972951251,
  32.247973266361804, 33.8561244442423, 35.46538285144842, 37.07554696139253, 38.68484536164441,
  40.29666519026153, 44.33041377334919,
};

// this is log10(partition function)
double Cr48_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 8.685880952436747e-07, 8.251516766996927e-06, 3.951899976600419e-05, 0.00013330794422173613,
  0.00035120219371925, 0.006401856055765157, 0.02685304570895992, 0.06215311825135839, 0.10696594975266842,
  0.15598699109465686, 0.205815844445829, 0.25471214514215257, 0.30198352738731143, 0.39152612205819926,
  0.47640596203905256, 0.5602400543128645, 0.6474755901642433, 0.7433846322638775, 0.983175072037813,
  1.3096301674258988, 1.7067177823367587, 2.1398790864012365, 2.577491799837225, 3.012837224705172,
  3.437750562820388, 3.8549130223078554, 4.264817823009537, 4.666517980554881, 5.648360010980932,
  6.606381365110605, 7.550228353055094, 8.484299839346786, 9.414973347970818, 10.340444114840118,
  11.264817823009537, 12.1846914308176, 13.103803720955957, 14.021189299069938, 14.935003151453655,
  15.846337112129806, 16.75511226639507, 17.66181268553726, 18.565847818673518, 19.468347330412158,
  20.369215857410143, 21.267171728403014, 22.161368002234976, 23.056904851336473, 23.94939000664491,
  24.840733234611807, 25.73078227566639, 26.619093330626743, 27.50785587169583, 28.394451680826215,
  29.281033367247726, 30.164352855784436, 31.049218022670182, 31.934498451243567, 33.70156798505593,
  35.46834733041216, 37.23299611039215, 38.99913054128737, 40.764922984649886, 42.5327543789925,
  44.30102999566398, 48.727541257028555,
};

// this is log10(partition function)
double Fe52_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 1.7371744532199383e-06, 9.554373504133797e-06, 3.778197643341552e-05,
  0.00011333607006293108, 0.0030242952161453874, 0.015422212189991184, 0.040215337130588114, 0.07478865660777631,
  0.11488541698288196, 0.15714990338033966, 0.19960737134331175, 0.24132628928072955, 0.3217032118192907,
  0.3993396534463543, 0.4778337814344742, 0.5623989859221217, 0.6594581913549248, 0.9153998352122699,
  1.2695129442179163, 1.6910814921229684, 2.143014800254095, 2.6009728956867484, 3.0569048513364727,
  3.503790683057181, 3.946452265013073, 4.383815365980431, 4.818225893613955, 5.888740960682893,
  6.944482672150168, 7.989894563718773, 9.02938377768521, 10.060697840353612, 11.086359830674748,
  12.11058971029925, 13.127104798364808, 14.139879086401237, 15.14921911265538, 16.152288344383056,
  17.152288344383056, 18.14921911265538, 19.143014800254097, 20.133538908370216, 21.12057393120585,
  22.103803720955955, 23.08635983067475, 24.06445798922692, 25.041392685158225, 26.01703333929878,
  26.989449817666692, 27.960946195733833, 28.930949031167522, 29.899273187317604, 30.8668778143375,
  31.833784374656478, 32.79934054945358, 33.76417613239033, 34.72835378202123, 36.655138434811384,
  38.58092497567562, 40.505149978319906, 42.42975228000241, 44.3541084391474, 46.28103336724773,
  48.20682587603185, 53.02938377768521,
};

// this is log10(partition function)
double Ni56_pf_array[] = {
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 4.342942647204277e-07, 7.817230319428646e-06, 6.42708273977769e-05,
  0.0002904458650804842, 0.0009123622824012837, 0.0022498876258026487, 0.0046944487518873, 0.014735532704563181,
  0.03529042138996706, 0.07190703372466718, 0.13162956968664008, 0.22190042758492473, 0.5092025223311029,
  0.9132839017604184, 1.374748346010104, 1.8555191556678001, 2.3404441148401185, 2.8221680793680175,
  3.3031960574204886, 3.783903579272735, 4.26245108973043, 4.7419390777291985, 5.9344984512435675,
  7.117271295655764, 8.292256071356476, 9.456366033129044, 10.608526033577194, 11.750508394851346,
  12.88309335857569, 14.008600171761918, 15.123851640967086, 16.232996110392154, 17.33645973384853,
  18.432969290874407, 19.525044807036846, 20.612783856719737, 21.695481676490196, 22.773786444981194,
  23.8481891169914, 24.919078092376076, 25.987219229908003, 27.053078443483418, 28.113943352306837,
  29.17609125905568, 30.232996110392154, 31.287801729930226, 32.3424226808222, 33.39619934709574,
  34.44715803134222, 35.49692964807321, 36.54530711646582, 37.594392550375424, 39.68752896121463,
  41.77959649125783, 43.86981820797933, 45.959518376973, 48.04921802267018, 50.13987908640124,
  52.230448921378276, 57.462397997898954,
};

void interpolate_pf( double t9 , double * pf_array , double * pf , double * dpf_dT ){

  if( t9 >= temp_array[0] && t9 < temp_array[npts_1-1] ){

    // find the largest temperature element <= t9 using a binary search

    int left = 0;
    int right = npts_1;

    while( left < right ){
      int mid = (left + right) / 2;
        if( temp_array[mid] > t9 ){
          right = mid;
        } else {
          left = mid + 1;
        }
    }

    int idx = right - 1;

    // now we have temp_array[idx] <= t9 < temp_array[idx+1]

    // construct the slope -- this is (log10(pf_{i+1}) - log10(pf_i)) / (T_{i+1} - T_i)

    double slope = (pf_array[idx+1] - pf_array[idx]) / (temp_array[idx+1] - temp_array[idx]);

    // find the PF

    double log10_pf = pf_array[idx] + slope * (t9 - temp_array[idx]);
    *pf = pow(10.0, log10_pf);

    double dpf_dT9 = pow(10.0, log10_pf) * M_LN10 * slope;
    *dpf_dT = dpf_dT9 * 1.e-9;

  } else {
    
    // T < the smallest T or >= the largest T in the partition function table
    *pf = 1.0;
    *dpf_dT = 0.0;
  }
}

void get_partition_function( int inuc , double T9 , double * pf , double * dpf_dT ){

  // inuc is the 1-based index for the species

  switch (inuc) {

  case O16:
    interpolate_pf( T9 , O16_pf_array , pf , dpf_dT );
    break;

  case Ne20:
    interpolate_pf( T9 , Ne20_pf_array , pf , dpf_dT );
    break;

  case Mg24:
    interpolate_pf( T9 , Mg24_pf_array , pf , dpf_dT );
    break;

  case Si28:
    interpolate_pf( T9 , Si28_pf_array , pf , dpf_dT );
    break;

  case S32:
    interpolate_pf( T9 , S32_pf_array , pf , dpf_dT );
    break;

  case Ar36:
    interpolate_pf( T9 , Ar36_pf_array , pf , dpf_dT );
    break;

  case Ca40:
    interpolate_pf( T9 , Ca40_pf_array , pf , dpf_dT );
    break;

  case Ti44:
    interpolate_pf( T9 , Ti44_pf_array , pf , dpf_dT );
    break;

  case Cr48:
    interpolate_pf( T9 , Cr48_pf_array , pf , dpf_dT );
    break;

  case Fe52:
    interpolate_pf( T9 , Fe52_pf_array , pf , dpf_dT );
    break;

  case Ni56:
    interpolate_pf( T9 , Ni56_pf_array , pf , dpf_dT );
    break;

  default:
    *pf = 1.0;
    *dpf_dT = 0.0;
  }
}

double get_spin_state( int inuc ){

    double spin = -1.0;

    switch (inuc) {

    case He4:
    case C12:
    case O16:
    case Ne20:
    case Mg24:
    case Si28:
    case S32:
    case Ar36:
    case Ca40:
    case Ti44:
    case Cr48:
    case Fe52:
    case Ni56:
        spin = 1;
        break;
    }

    return spin;

}

//////////////////////////////////////////////////////////////////////////////////

void get_T_independent_nse_state( struct NSE * state ){
  // Computes (2J + 1) m^{5/2} (kB 1e9 / 2 π ħ^2)^{3/2} / ρ
  // a portion of the nse_state that is not dependent on Temperature.

  double kBGK_2PiHbar2 = (k_B * 1.0e9) / (2.0 * M_PI * hbar * hbar);
  double rho_inv = 1.0 / state->rho;
  double spin = 1.0;

  for( int n=0 ; n<NUM_I ; ++n ){
    // Absorb mion into kBGK_2PiHbar2 to prevent underflow at compile time.
    spin = get_spin_state(n);
    state->xn[n] = mion[n] * spin * pow(mion[n]*kBGK_2PiHbar2, 1.5) * rho_inv;
  }
}

void compute_coulomb_contribution( double u_c[] , double du_c_dT9[] , struct NSE * state ){
  // This function computes the coulomb contribution and its T9 derivative

  // if we use chabrier1998 screening
  // Get the required terms to calculate coulomb correction term, u_c

  double T_in = state->T;

  // Find n_e for original state;
  // Note that y_e depends on the mass fraction,
  // but we use the coulomb correction to compute the mass fraction
  // But the one of the constraint is on input y_e
  // So here y_e is simply the actual y_e we want to achieve.
  // so u_c only depends on temperature.

  double T9i = 1.0e9 / T_in;
  double n_e = state->rho * state->y_e * n_A;
  double gamma_e_constants = q_e * q_e * pow(4.0 * M_PI / 3.0, 1.0 / 3.0) / (k_B * 1.0e9);
  double Gamma_e = gamma_e_constants * cbrt(n_e) * T9i;

  for( int n=0 ; n<NUM_I ; ++n ){
    // term for calculating u_c
    double Z_53 = pow(zion[n], 5.0 / 3.0);
    double Gamma = Z_53 * Gamma_e;

    // chemical potential for coulomb correction
    // see appendix of Calder 2007, doi:10.1086/510709 for more detail

    // reuse existing implementation from screening routine

    //
    // Here u_c is a dimensionless quantity.
    // i.e. normalized by kB T
    //

    // Fitted parameters, see Chabrier & Potekhin 1998 Sec.IV
    double A_1 = -0.9052;
    double A_2 = 0.6322;
    double A_3 = -0.5 * sqrt(3.0) - A_1 / sqrt(A_2);

    // Precompute some expressions that are reused in the derivative
    double sqrt_Gamma = sqrt(Gamma);
    double sqrt_1_Gamma_A2 = sqrt(1.0 + Gamma/A_2);
    double sqrt_Gamma_A2_Gamma = sqrt(Gamma * (A_2 + Gamma));
    double sqrt_Gamma_A2 = sqrt(Gamma/A_2);

    double f = A_1 * (sqrt_Gamma_A2_Gamma - A_2 * log(sqrt_Gamma_A2 + sqrt_1_Gamma_A2)) + 2.0 * A_3 * (sqrt_Gamma - atan(sqrt_Gamma));
    //u_c[n] = k_B * T_in / (MeV2eV * ev2erg) * f;
    u_c[n] = f;

    // Find dΓ / dT9 = - Γ / T9
    double dGamma_dT9 = - Gamma * T9i;

    // Find temperature derivative of the coulomb potential term: d (μ^c/kBT) / dT
    du_c_dT9[n] = dGamma_dT9 * sqrt(Gamma) * (A_3 / (1.0 + Gamma) + A_1 / sqrt(A_2 + Gamma));
  }
}

void get_T_dependent_nse_state( struct NSE * state ){
  // It does two things:
  // 1. Applies T9^{3/2} to the input state,
  //    where the input state should come from get_T_independent_nse_state()
  // 2. Evaluates the temperature dependent piece in the exponent
  //    i.e. - μ^c/kBT  + ln(pf)
  //    This is the term will be used inside the exp()
  //    It also evaluates the dlog/T9 derivative, including the effect of T9^{3/2}
  //    i.e. dlog/dT9 = 1.5 / T9 - du^c/dT9 + dlogpf/dT9

  double T_in = state->T;
  double T9 = T_in * 1.0e-9;
  double T9_32 = sqrt(T9*T9*T9);
  double T9i = 1.0 / T9;

  // Evaluate the coulomb potential: μ^c / kBT and potentially find its T9 derivative
  double u_c[NUM_I];
  double du_c_dT9[NUM_I];
  //if( SCREEN_METHOD == SCREEN_METHOD_chabrier1998 ){
  //compute_coulomb_contribution( u_c , du_c_dT9 , state );
  for( int n=0 ; n<NUM_I ; ++n ){ u_c[n]=0; du_c_dT9[n]=0; }
  //}

  for( int n=0 ; n<NUM_I ; ++n ){

    // Evaluate internal nuclear partition function, in log space

    double pf;
    double dlogpf_dT9;
    get_partition_function( n , T9 , &pf , &dlogpf_dT9 );

    // Compute the log term
    state->log_data[n] = log(pf) - u_c[n];
    // Apply T9^{3/2} to the state
    state->xn[n] *= T9_32;

    state->dlogX_dT9[n] = 1.5 * T9i + dlogpf_dT9 - du_c_dT9[n];
  }
}

void apply_nse_exponent( struct NSE * state ){
    // Last step to compute the NSE mass fractions.
    // It applies the nse exponent which depends on
    // the chemical potential and the precomputed Temperature term
    // from get_T_dependent_exponent_nse_state

    // This also completes dlogX_dT9 with the contribution from:
    // (Z_k μ_p + N_k μ_n + B_k) / (kB T)

    double T_in = state->T;
    double T9i = 1.0e9 / T_in;
    double ikMeV_GK = 1.0e-9 * (MeV2eV * ev2erg) / k_B;
    double ikTMeV = ikMeV_GK * T9i;

    for( int n=0 ; n<NUM_I ; ++n ){

        // Apply the exponent part

        // prevent an overflow on exp by capping the exponent -- we hope that a subsequent
        // iteration will make it happy again
        double mu_kBT = (zion[n] * state->mu_p + (aion[n] - zion[n]) * state->mu_n + bion[n]) * ikTMeV;
        double exponent = fmin(500.0, state->log_data[n] + mu_kBT);
        state->xn_exp[n] = state->xn[n]*exp(exponent);

        state->dlogX_dT9[n] += -mu_kBT * T9i;
    }

    // Compute the corresponding electron fraction, Ye.

    state->y_e = 0.0;
    for( int n=0 ; n<NUM_I ; ++n ){
        state->y_e += state->xn_exp[n] * zion[n] * aion_inv[n];
    }
}

////////////////////////////////////////////////////////////////////////////////////

// constraint equation
void fcn( double x[] , double fvec[] , struct NSE * state ){
  
  // Apply exponent component for calculating nse mass fractions
  state->mu_p = x[0];
  state->mu_n = x[1];
  apply_nse_exponent( state );

  // constraint equation 1, mass fraction sum to 1
  fvec[0] = -1.0;
  for( int n=0 ; n<NUM_I ; ++n ){
    fvec[0] += state->xn_exp[n];
  }

  // constraint equation 2, electron fraction should be the same
  fvec[1] = state->y_e - state->y_e_init;

}

// constraint jacobian
void jcn( double x[] , double fjac[][2] , struct NSE * state ){
    
  // Apply exponent component for calculating nse mass fractions
  state->mu_p = x[0];
  state->mu_n = x[1];
  apply_nse_exponent( state );

    // evaluate jacobian of the constraint
    fjac[0][0] = 0.0;
    fjac[0][1] = 0.0;
    fjac[1][0] = 0.0;
    fjac[1][1] = 0.0;

    double T_in = state->T;

    for( int n=0 ; n<NUM_I ; ++n ){
        fjac[0][0] += state->xn_exp[n] * zion[n] / k_B / T_in * (MeV2eV * ev2erg);
        fjac[0][1] += state->xn_exp[n] * (aion[n] - zion[n]) / k_B / T_in * (MeV2eV * ev2erg);
        fjac[1][0] += state->xn_exp[n] * zion[n] * zion[n] * aion_inv[n] / k_B / T_in * (MeV2eV * ev2erg);
        fjac[1][1] += state->xn_exp[n] * zion[n] * (aion[n] - zion[n]) * aion_inv[n] / k_B / T_in * (MeV2eV * ev2erg);
    }
}

// the six hybrj solver functions

double enorm( int n , double * x ){
  // int n = neqs (number of equations)
  // double * x = f[neqs] (fvec in Castro terms)

  // given an n-vector x, this function calculates the
  // euclidean norm of x.
  //
  // the euclidean norm is computed by accumulating the sum of
  // squares in three different sums. the sums of squares for the
  // small and large components are scaled so that no overflows
  // occur. non-destructive underflows are permitted. underflows
  // and overflows do not occur in the computation of the unscaled
  // sum of squares for the intermediate components.
  // the definitions of small, intermediate and large components
  // depend on two constants, rdwarf and rgiant. the main
  // restrictions on these constants are that rdwarf**2 not
  // underflow and rgiant**2 not overflow. the constants
  // given here are suitable for every known computer.
  //
  // argonne national laboratory. minpack project. march 1980.
  // burton s. garbow, kenneth e. hillstrom, jorge j. more

  double rdwarf = 3.834e-20;
  double rgiant = 1.304e19;

  double s1 = 0.0;
  double s2 = 0.0;
  double s3 = 0.0;

  double x1max = 0.0;
  double x3max = 0.0;

  double agiant = rgiant / (double)n;

  for( int i=0 ; i<n ; ++i ){
    double xabs = fabs(x[i]);

    if( isnan(xabs) || isinf(xabs) ){
      return 0.0;
    }
    if( xabs <= rdwarf || xabs >= agiant ){
      if( xabs > rdwarf ){

        // sum for large components.

        if( xabs > x1max ){
          s1 = 1.0 + s1 * pow(x1max/xabs,2.0);
          x1max = xabs;
        } else {
          s1 += pow(xabs/x1max,2.0);
        }
      } else {

        // sum for small components.

        if( xabs > x3max ){
          s3 = 1.0 + s3 * pow(x3max/xabs,2.0);
          x3max = xabs;
        } else {
          if (xabs != 0.0) {
            s3 += pow(xabs/x3max,2.0);
          }
        }
      }
    } else {

      // sum for intermediate components.

      s2 += xabs * xabs;
    }
  }

  // calculation of norm.

  double _enorm = 0.0;
  if( s1 != 0.0 ){
    _enorm = x1max * sqrt(s1 + (s2/x1max) / x1max);
  } else {
    if( s2 != 0.0 ){
      if( s2 >= x3max ){
        _enorm = sqrt(s2 * (1.0 + (x3max/s2)*(x3max*s3)));
      }
      if( s2 < x3max ){
        _enorm = sqrt(x3max * ((s2/x3max) + (x3max*s3)));
      }
    } else {
      _enorm = x3max * sqrt(s3);
    }
  }

  return _enorm;
}

void qrfac( int neqs , double a[][neqs] , double * rdiag , double * acnorm , double * wa ){

  // temporary storage for a column of the Jacobian -- the Fortran
  // code would pass to enorm something like a(i, j) to compute the
  // norm of column j from row i to the last row.

  // this subroutine uses householder transformations with column
  // pivoting (optional) to compute a qr factorization of the
  // m by n matrix a. that is, qrfac determines an orthogonal
  // matrix q, a permutation matrix p, and an upper trapezoidal
  // matrix r with diagonal elements of nonincreasing magnitude,
  // such that a*p = q*r. the householder transformation for
  // column k, k = 1,2,...,min(m,n), is of the form
  //                        t
  //        i - (1/u(k))*u*u
  //
  // where u has zeros in the first k-1 positions. the form of
  // this transformation and the method of pivoting first
  // appeared in the corresponding linpack subroutine.

  // changes from the Fortran version:
  //
  // * removed the ability to pivot
  // * assume input is square
  // * written to use our hybrj_t type

  // jac is an n by n array. on input a contains the matrix for
  // which the qr factorization is to be computed. on output
  // the strict upper trapezoidal part of a contains the strict
  // upper trapezoidal part of r, and the lower trapezoidal
  // part of a contains a factored form of q (the non-trivial
  // elements of the u vectors described above).
  //
  // the following are assigned to different work arrays in the hybrj type
  //
  // rdiag is an output array of length n which contains the
  // diagonal elements of r.
  //
  // acnorm is an output array of length n which contains the
  // norms of the corresponding columns of the input matrix a.
  // if this information is not needed, then acnorm can coincide
  // with rdiag.
  //
  // wa is a work array of length n. if pivot is false, then wa
  // can coincide with rdiag.
  //
  //
  // argonne national laboratory. minpack project. march 1980.
  // burton s. garbow, kenneth e. hillstrom, jorge j. more

  double tmp[neqs];

  // compute the initial column norms and initialize several arrays.

  for( int j=0 ; j<neqs ; ++j ){
    for( int irow=0 ; irow<neqs ; ++irow ){
      tmp[irow] = a[irow][j];
    }
    acnorm[j] = enorm(neqs, tmp);
    rdiag[j] = acnorm[j];
    wa[j] = rdiag[j];
  }

  // reduce a to r with householder transformations.

  for( int j=0 ; j<neqs ; ++j ){

    // compute the householder transformation to reduce the
    // j-th column of a to a multiple of the j-th unit vector.

    // we want the norm of a(j,j:n)
    for( int irow=j ; irow<neqs ; ++irow ){
      tmp[irow-j] = a[irow][j];
    }
    double ajnorm = enorm(neqs-j, tmp);
    if( ajnorm != 0.0 ){
      if( a[j][j] < 0.0 ){
        ajnorm = -ajnorm;
      }
      for( int irow=j ; irow<neqs ; ++irow ){
        a[irow][j] /= ajnorm;
      }
      a[j][j] += 1.0;

      // apply the transformation to the remaining columns
      // and update the norms.

      int jp1 = j + 1;
      if( neqs > jp1 ){
        for( int k=jp1 ; k<neqs ; ++k ){
          double sum = 0.0;
          for( int i=j ; i<neqs ; ++i ){
            sum += a[i][j] * a[i][k];
          }
          double temp = sum / a[j][j];
          for( int i=j ; i<neqs ; ++i ){
            a[i][k] += -temp * a[i][j];
          }
        }
      }
    }
    rdiag[j] = -ajnorm;
  }

}

void qform( int neqs , double q[][neqs], double * wa ){
  // q is an n by n array. on input the full lower trapezoid in
  // the n columns of q contains the factored form.
  // on output q has been accumulated into a square matrix.

  // this subroutine proceeds from the computed qr factorization of
  // an m by n matrix a to accumulate the m by m orthogonal matrix
  // q from its factored form.
  //
  //  argonne national laboratory. minpack project. march 1980.
  //  burton s. garbow, kenneth e. hillstrom, jorge j. more

  // changes:
  // originally, this routine supported an mxn matrix, but we assume m=n

  // zero out upper triangle of q in the first n

  if( neqs >= 2 ){
    for( int j=1 ; j<neqs ; ++j ){
      int jm1 = j - 1;
      for( int i=0 ; i<=jm1 ; ++i ){
        q[i][j] = 0.0;
      }
    }
  }

  // accumulate q from its factored form.

  for( int l=0 ; l<neqs ; ++l ){
    int k = neqs - l - 1;
    for( int i=k ; i<neqs ; ++i ){
      wa[i] = q[i][k];
      q[i][k] = 0.0;
    }
    q[k][k] = 1.0;
    if( wa[k] != 0.0 ){
      for( int j=k ; j<neqs ; ++j ){
        double sum = 0.0;
        for( int i=k ; i<neqs ; ++i ){
          sum += q[i][j] * wa[i];
        }
        double temp = sum / wa[k];
        for( int i=k ; i<neqs ; ++i ){
          q[i][j] += -temp * wa[i];
        }
      }
    }
  }
}

void dogleg( int neqs , double * r , double * diag , double * qtb , double delta , double * x , double * wa1 , double * wa2 ){
            
  // given an n by n matrix a, an n by n nonsingular diagonal
  // matrix d, an n-vector b, and a positive number delta, the
  // problem is to determine the convex combination x of the
  // gauss-newton and scaled gradient directions that minimizes
  // (a*x - b) in the least squares sense, subject to the
  // restriction that the euclidean norm of d*x be at most delta.
  //
  // this subroutine completes the solution of the problem
  // if it is provided with the necessary information from the
  // qr factorization of a. that is, if a = q*r, where q has
  // orthogonal columns and r is an upper triangular matrix,
  // then dogleg expects the full upper triangle of r and
  // the first n components of (q transpose)*b.
  //
  // parameters:
  //
  //  r is an input array of length (n*(n+1))/2 which must contain the upper
  //    triangular matrix r stored by rows.
  //
  //  diag is an input array of length n which must contain the
  //    diagonal elements of the matrix d.
  //
  //  qtb is an input array of length n which must contain the first
  //    n elements of the vector (q transpose)*b.
  //
  //  delta is a positive input variable which specifies an upper
  //    bound on the euclidean norm of d*x.
  //
  //  x is an output array of length n which contains the desired
  //    convex combination of the gauss-newton direction and the
  //    scaled gradient direction.
  //
  //  wa1 and wa2 are work arrays of length n.
  //
  //  argonne national laboratory. minpack project. march 1980.
  //  burton s. garbow, kenneth e. hillstrom, jorge j. more

  double epsmch = DBL_EPSILON;

  // first, calculate the gauss-newton direction.

  int jj = (neqs*(neqs + 1)) / 2;
  for( int k=0 ; k<neqs ; ++k ){
    int j = neqs - k - 1;
    int jp1 = j + 1;
    jj = jj - k - 1;
    int l = jj + 1;

    double sum = 0.0;
    if( neqs > jp1 ){
      for( int i = jp1 ; i<neqs ; ++i ){
        sum += r[l] * x[i];
        l++;
      }
    }
    double temp = r[jj];
    if( temp == 0.0 ){
      l = j;
      for( int i=0 ; i<=j ; ++i ){
        temp = fmax(temp, fabs(r[l]));
        l = l + neqs - 1 - i;
      }
      temp *= epsmch;
      if( temp == 0.0 ){
        temp = epsmch;
      }
    }

    x[j] = (qtb[j] - sum) / temp;
  }

  // test whether the gauss-newton direction is acceptable.

  for( int j=0 ; j<neqs ; ++j ){
    wa1[j] = 0.0;
    wa2[j] = diag[j] * x[j];
  }
  double qnorm = enorm(neqs, wa2);

  if( qnorm <= delta ){ return; }

  // the gauss-newton direction is not acceptable.
  // next, calculate the scaled gradient direction.

  int l = 0;
  for( int j=0 ; j<neqs ; ++j ){
    double temp = qtb[j];
    for( int i=j ; i<neqs ; ++i ){
      wa1[i] += r[l] * temp;
      l++;
      if( isnan(wa1[i]) || isinf(wa1[i]) ){
        return;
      }
    }
    wa1[j] /= diag[j];

    if( isnan(wa1[j]) || isinf(wa1[j]) ){
      return;
    }
  }

  // calculate the norm of the scaled gradient and test for
  // the special case in which the scaled gradient is zero.

  double gnorm = enorm(neqs, wa1);
  double sgnorm = 0.0;
  double alpha = delta / qnorm;

  if( gnorm != 0.0 ){

    // calculate the point along the scaled gradient
    // at which the quadratic is minimized.

    for( int j=0 ; j<neqs ; ++j ){
      wa1[j] = (wa1[j] / gnorm) / diag[j];
    }

    l = 0;
    for( int j=0 ; j<neqs ; ++j ){
      double sum = 0.0;
      for( int i=j ; i<neqs ; ++i ){
        sum += r[l] * wa1[i];
        l++;
      }
      wa2[j] = sum;
    }
    double temp = enorm(neqs, wa2);
    sgnorm = (gnorm / temp) / temp;

    // test whether the scaled gradient direction is acceptable.

    alpha = 0.0;
    if( sgnorm < delta ){

      // the scaled gradient direction is not acceptable.
      // finally, calculate the point along the dogleg
      // at which the quadratic is minimized.

      double bnorm = enorm(neqs, qtb);
      temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta);
      temp = temp - (delta/qnorm) * pow(sgnorm/delta,2.0) + sqrt( pow(temp-(delta/qnorm),2.0) + (1.0-pow(delta/qnorm,2.0)) * (1.0-pow(sgnorm/delta,2.0)));
      alpha = ((delta/qnorm)*(1.0 - pow(sgnorm/delta,2.0))) / temp;
    }
  }

  // form appropriate convex combination of the gauss-newton
  // direction and the scaled gradient direction.

  double temp = (1.0 - alpha) * fmin(sgnorm, delta);

  for( int j=0 ; j<neqs ; ++j ){
    x[j] = temp * wa1[j] + alpha * x[j];
  }

}

void r1updt( int neqs , double * s , double * u , double * v , double * w , int * sing ){
           
  // given an n by n lower trapezoidal matrix s, an n-vector u,
  // and an n-vector v, the problem is to determine an
  // orthogonal matrix q such that
  //
  //                   t
  //           (s + u*v )*q
  //
  // is again lower trapezoidal.
  //
  // this subroutine determines q as the product of 2*(n - 1)
  // transformations
  //
  //          gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
  //
  // where gv(i), gw(i) are givens rotations in the (i,n) plane
  // which eliminate elements in the i-th and n-th planes,
  // respectively. q itself is not accumulated, rather the
  // information to recover the gv, gw rotations is returned.

  // changes:
  // we now assume that the matrix is nxn
  //
  // s is an array of length (n*(n+1))/2. on input s must contain the lower
  //   trapezoidal matrix s stored by columns. on output s contains
  //   the lower trapezoidal matrix produced as described above.
  //
  // u is an input array of length m which must contain the
  //   vector u.
  //
  // v is an array of length n. on input v must contain the vector
  //   v. on output v(i) contains the information necessary to
  //   recover the givens rotation gv(i) described above.
  //
  //  w is an output array of length m. w(i) contains information
  //    necessary to recover the givens rotation gw(i) described
  //    above.
  //
  //  sing is a logical output variable. sing is set true if any
  //       of the diagonal elements of the output s are zero. otherwise
  //       sing is set false.
  //
  // argonne national laboratory. minpack project. march 1980.
  // burton s. garbow, kenneth e. hillstrom, jorge j. more,
  // john l. nazareth

  double giant = DBL_MAX;

  // initialize the diagonal element pointer.

  int jj = (neqs * (neqs + 1))/2;

  // move the nontrivial part of the last column of s into w.

  int l = jj;
  for( int i=neqs ; i <= neqs; ++i) {
    w[i-1] = s[l-1];
    l++;
  }

  // rotate the vector v into a multiple of the n-th unit vector
  // in such a way that a spike is introduced into w.

  int nm1 = neqs - 1;
  if( nm1 >= 1 ){
    for( int nmj=1 ; nmj<=nm1 ; ++nmj ){
      int j = neqs - 1 - nmj;
      jj += -(neqs - j);
      w[j] = 0.0;
      if( v[j] != 0.0 ){

        // determine a givens rotation which eliminates the
        // j-th element of v.

        double tau;
        double fcos;
        double fsin;

        if( fabs(v[neqs-1]) < fabs(v[j]) ){
          double cotan = v[neqs-1] / v[j];
          fsin = 0.5 / sqrt(0.25 + 0.25 * cotan * cotan);
          fcos = fsin * cotan;
          tau = 1.0;
          if( fabs(fcos) * giant > 1.0 ){
            tau = 1.0 / fcos;
          }
        } else {
          double ftan = v[j] / v[neqs-1];
          fcos = 0.5 / sqrt(0.25 + 0.25 * ftan * ftan);
          fsin = fcos * ftan;
          tau = fsin;
        }

        // apply the transformation to v and store the information
        // necessary to recover the givens rotation.

        v[neqs-1] = fsin * v[j] + fcos * v[neqs-1];
        v[j] = tau;

        // apply the transformation to s and extend the spike in w.

        l = jj;
        for( int i=j ; i<neqs ; ++i ){
          double temp = fcos * s[l] - fsin * w[i];
          w[i] = fsin * s[l] + fcos * w[i];
          s[l] = temp;
          l++;
        }
      }
    }
  }

  // add the spike from the rank 1 update to w.

  for( int i=0 ; i<neqs ; ++i ){
    w[i] += v[neqs-1] * u[i];
  }

  // eliminate the spike.

  *sing = 0; //false

  if( nm1 >= 1 ){
    for( int j=0 ; j<nm1 ; ++j ){
      if( w[j] != 0.0 ){

        // determine a givens rotation which eliminates the
        // j-th element of the spike.

        double tau;
        double fcos;
        double fsin;

        if( fabs(s[jj]) < fabs(w[j]) ){
          double cotan = s[jj] / w[j];
          fsin = 0.5 / sqrt(0.25 + 0.25 * cotan * cotan);
          fcos = fsin * cotan;
          tau = 1.0;
          if( fabs(fcos) * giant > 1.0 ){
            tau = 1.0 / fcos;
          }
        } else {
          double ftan = w[j] / s[jj];
          fcos = 0.5 / sqrt(0.25 + 0.25 * ftan * ftan);
          fsin = fcos * ftan;
          tau = fsin;
        }

        // apply the transformation to s and reduce the spike in w.

        l = jj;
        for( int i=j ; i<neqs ; ++i ){
          double temp = fcos * s[l] + fsin * w[i];
          w[i] = -fsin * s[l] + fcos * w[i];
          s[l] = temp;
          l++;
        }

        // store the information necessary to recover the
        // givens rotation.

        w[j] = tau;
      }

      // test for zero diagonal elements in the output s.

      if( s[jj] == 0.0 ){
        *sing = 1; //true
      }
      jj += (neqs - j);
    }
  }

  // move w back into the last column of the output s.

  l = jj;
  for( int i=neqs ; i<=neqs ; ++i ){
    s[l-1] = w[i-1];
    l++;
  }

  if( s[jj] == 0.0 ){
    *sing = 1; //true
  }
}

// given an m by n matrix a, this subroutine computes a*q where
//  q is the product of 2*(n - 1) transformations
//
//       gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
//
// and gv(i), gw(i) are givens rotations in the (i,n) plane which
// eliminate elements in the i-th and n-th planes, respectively.
// q itself is not given, rather the information to recover the
// gv, gw rotations is supplied.
//
// changes:
// * Split into 2 versions, one takes an nxn matrix and one takes a vector of length n

void r1mpyq_nxn( int neqs , double a[][neqs] , double * v , double * w ){

  // a is an neqs by neqs array. on input a must contain the matrix
  //   to be postmultiplied by the orthogonal matrix q
  //   described above. on output a*q has replaced a.
  //
  // v is an input array of length neqs. v(i) must contain the
  //   information necessary to recover the givens rotation gv(i)
  //   described above.
  //
  // w is an input array of length neqs. w(i) must contain the
  //   information necessary to recover the givens rotation gw(i)
  //   described above.
  //
  // argonne national laboratory. minpack project. march 1980.
  // burton s. garbow, kenneth e. hillstrom, jorge j. more

  // apply the first set of givens rotations to a.

  int nm1 = neqs - 1;

  double fcos;
  double fsin;

  if( nm1 >= 1 ){
    for( int nmj=2 ; nmj<=neqs ; ++nmj ){
      int j = neqs - nmj;
      if( fabs(v[j]) > 1.0 ){
        fcos = 1.0 / v[j];
        fsin = sqrt(1.0 - fcos*fcos);
      } else {
        // if( fabs(v[j]) <= 1.0_rt ){
        fsin = v[j];
        fcos = sqrt(1.0 - fsin*fsin);
      }
      for( int i=0 ; i<neqs ; ++i ){
        double temp = fcos * a[i][j] - fsin * a[i][nm1];
        a[i][nm1] = fsin * a[i][j] + fcos * a[i][nm1];
        a[i][j] = temp;
      }
    }

    // apply the second set of givens rotations to a.

    for( int j=0 ; j<nm1 ; ++j ){
      if( fabs(w[j]) > 1.0 ){
        fcos = 1.0 / w[j];
        fsin = sqrt(1.0 - fcos*fcos);
      } else {
      // if( fabs(w[j]) <= 1.0 ){
        fsin = w[j];
        fcos = sqrt(1.0 - fsin*fsin);
      }
      for( int i=0 ; i<neqs ; ++i ){
        double temp = fcos * a[i][j] + fsin * a[i][nm1];
        a[i][nm1] = -fsin*a[i][j] + fcos*a[i][nm1];
        a[i][j] = temp;
      }
    }
  }

}

void r1mpyq( int neqs , double * a , double * v , double * w ){

  // a is an neqs vector.  This is based on the original version
  // with m = 1.
  //
  // v is an input array of length neqs. v(i) must contain the
  //   information necessary to recover the givens rotation gv(i)
  //   described above.
  //
  // w is an input array of length neqs. w(i) must contain the
  //   information necessary to recover the givens rotation gw(i)
  //   described above.
  //
  // argonne national laboratory. minpack project. march 1980.
  // burton s. garbow, kenneth e. hillstrom, jorge j. more

  // apply the first set of givens rotations to a.


  int nm1 = neqs - 1;

  double fcos;
  double fsin;

  if( nm1 >= 1 ){
    for( int nmj=2 ; nmj<=neqs ; ++nmj ){
      int j = neqs - nmj;
      if( fabs(v[j]) > 1.0 ){
        fcos = 1.0 / v[j];
        fsin = sqrt(1.0 - fcos*fcos);
      } else {
        // if( fabs(v[j]) <= 1.0 ){
        fsin = v[j];
        fcos = sqrt(1.0 - fsin*fsin);
      }
      double temp = fcos*a[j] - fsin*a[nm1];
      a[nm1] = fsin*a[j] + fcos*a[nm1];
      a[j] = temp;
    }

    // apply the second set of givens rotations to a.

    for( int j=0 ; j<nm1 ; ++j ){
      if( fabs(w[j]) > 1.0 ){
        fcos = 1.0 / w[j];
        fsin = sqrt(1.0 - fcos*fcos);
      } else {
        // if( fabs(w[j]) <= 1.0 ){
        fsin = w[j];
        fcos = sqrt(1.0 - fsin*fsin);
      }
      double temp = fcos*a[j] + fsin*a[nm1];
      a[nm1] = -fsin*a[j] + fcos*a[nm1];
      a[j] = temp;
    }
  }
}

// A Powell-Hybrid solver for finding nse state used for calibrating
// chemical potential of proton and neutron
void nse_hybrid_solver( int neqs , struct NSE * state , double eps ){

  // Set relative error between two iterations
  double xtol = eps;

  // Different modes for the hybrid solver
  // mode = 1: scales internally
  // mode = 2: scales based on user-supplied diag

  // Choose mode 1 for all.
  // If T is unknown, we manually scale to work with T9
  int mode = 1;
  int info;
  double factor = 100.0;
  int maxfev = 1000;

  // random flag number for evaluation in for loop;

  int flag = 0;

  // Fine-tune variables

  double dx;
  int is_pos_new;
  int is_pos_old = 0; //false

  double f[neqs];
  double jac[neqs][neqs];
  int nfev = 0;
  int njev = 0;
  double x[neqs];
  double outer_x[neqs];
  double inner_x[neqs];
  double diag[neqs];
  double qtf[neqs];
  double r[neqs*(neqs+1)/2];
  double wa1[neqs];
  double wa2[neqs];
  double wa3[neqs];
  double wa4[neqs];

  outer_x[0] = state->mu_p;
  outer_x[1] = state->mu_n;
  if( neqs == 3 ){
    // For temperature, solve for T9.
    outer_x[2] = state->T * 1.e-9;
  }

  // fine tuning initial guesses
  // Here we only fine tune chemical potentials

  for( int i=0 ; i<20 ; ++i ){

    dx = 0.5;
    for( int n=0 ; n<neqs ; ++n ){
      inner_x[n] = outer_x[n];
    }

    for( int j=0 ; j<20 ; ++j ){

      for( int n=0 ; n<neqs ; ++n ){
        x[n] = inner_x[n];
      }

////////////////////////////////////////////////////////////////////////////////////
      //hybrj<neqs, nse_solver_data<T>>(hj, state_data,nse_fcn<nse_solver_data<T>>,nse_jcn<nse_solver_data<T>>);

      // the purpose of hybrj is to find a zero of a system of
      // n nonlinear functions in n variables by a modification
      // of the powell hybrid method. the user must provide a
      // subroutine which calculates the functions and the jacobian.
      //
      // fcn jcn are the names of the user-supplied subroutine which
      // calculates the functions and the jacobian.
      //
      // void fcn(Array1D<Real, 1, neqs>& x, Array1D<Real, 1, neqs>& fvec, int& iflag)
      // void jcn(Array1D<Real, 1, neqs>& x, Array2D<Real, 1, neqs, 1, neqs> fjac, int& iflag)
      //
      // the value of iflag should not be changed by fcn unless
      // the user wants to terminate execution of hybrj.
      // in this case set iflag to a negative integer.

      int finished = 0; //false

      double epsmch = DBL_EPSILON;

      info = 0;
      int iflag = 0;

      double delta = 0.0;
      double xnorm = 0.0;

      // check the input parameters for errors.
      if( xtol < 0.0 ){
        printf("error: xtol must be > 0");
      }

      if( mode == 2 ){
        for( int j=0 ; j<neqs ; ++j ){
          if( diag[j] <= 0.0 ){
            finished = 1;
            break;
          }
        }
      }

      // evaluate the function at the starting point
      // and calculate its norm.
      double fnorm = 0.0;

      iflag = 1;

      fcn(x, f, state);
      nfev = 1;
      if( iflag < 0 ){
        // user requested abort
        finished = 1;
      } else {
        fnorm = enorm(neqs, f);
      }

      // initialize iteration counter and monitors.

      int iter = 1;
      int ncsuc = 0;
      int ncfail = 0;
      int nslow1 = 0;
      int nslow2 = 0;

      // beginning of the outer loop.

      while( finished == 0 ){

        int jeval = 1; //true

        // calculate the jacobian matrix.

        jcn(x, jac, state);
        njev++;
        if( iflag < 0 ){
          finished = 1;
          break;
        }

        // compute the qr factorization of the jacobian.

        qrfac(neqs, jac, wa1, wa2, wa3);

        // on the first iteration, 
        // scale according to the norms of the columns of the initial jacobian,
        // calculate the norm of the scaled x,
        // and initialize the step bound delta.

        if( iter == 1 ){
          if( mode != 2 ){
            for( int j=0 ; j<neqs ; ++j ){
              diag[j] = wa2[j];
              if( wa2[j] == 0.0 ){
                diag[j] = 1.0;
              }
            }
          }

          for( int j=0 ; j<neqs ; ++j ){
            wa3[j] = diag[j] * x[j];
          }
          xnorm = enorm(neqs, wa3);
          delta = factor * xnorm;
          if( delta == 0.0 ){ delta = factor; }
        }

        // form (q transpose)*fvec and store in qtf.

        for( int i=0 ; i<neqs ; ++i ){
          qtf[i] = f[i];
        }
        for( int j=0 ; j<neqs ; ++j ){
          if( jac[j][j] != 0.0 ){
            double sum = 0.0;
            for( int i=j ; i<neqs ; ++i ){ 
              sum += jac[i][j] * qtf[i]; 
            }
            double temp = -sum / jac[j][j];
            for( int i=j ; i<neqs ; ++i ){
              qtf[i] += jac[i][j] * temp;
            }
          }
        }

        // copy the triangular factor of the qr factorization into r.

        int sing = 0; //false

        for( int j=0 ; j<neqs ; ++j ){
          int l = j;
          int jm1 = j-1;
          if( jm1 >= 0 ){
            for( int i=0; i<=jm1 ; ++i ){
              r[l] = jac[i][j];
              l = l + neqs - i;
            }
          }
          r[l] = wa1[j];
          if( wa1[j] == 0.0 ){ 
            sing = 1; 
          }
        }

        // accumulate the orthogonal factor in fjac.

        qform(neqs, jac, wa1);

        // rescale 
        if( mode != 2 ){
          for( int j=0 ; j<neqs ; ++j ){ 
            diag[j] = fmax( diag[j] , wa2[j] ); 
          }
        }

        // beginning of the inner loop.

        while( 1==1 ){

          // determine the direction p.

          dogleg(neqs, r, diag, qtf, delta, wa1, wa2, wa3);

          // store the direction p and x + p. calculate the norm of p.

          int invalid = 0; //false

          for( int i=0 ; i<neqs ; ++i ){
            if( isinf(wa1[i]) || isnan(wa1[i]) ){
              invalid = 1;
            }
          }

          if( invalid == 1 ){
            finished = 1;
            break;
          }

          for( int j=0 ; j<neqs ; ++j ){
            wa1[j] = -wa1[j];
            wa2[j] = x[j] + wa1[j];
            wa3[j] = diag[j] * wa1[j];
          }

          double pnorm = enorm(neqs, wa3);

          // on the first iteration, adjust the initial step bound.

          if( iter == 1 ){
            delta = fmin(delta, pnorm);
          }

          // evaluate the function at x + p and calculate its norm.

          fcn(wa2, wa4, state);
          nfev++;

          if( iflag < 0 ){
            finished = 1;
            break;
          }

          double fnorm1 = enorm(neqs, wa4);

          // compute the scaled actual reduction.

          double actred = -1.0;
          if( fnorm1 < fnorm ){
            actred = 1.0 - pow(fnorm1/fnorm,2.0);
          }

          // compute the scaled predicted reduction.

          int l = 0;
          for( int i=0 ; i<neqs ; ++i ){
            double sum = 0.0;
            for( int j=i ; j<neqs ; ++j ){
              sum += r[l] * wa1[j];
              l += 1;
            }
            wa3[i] = qtf[i] + sum;
          }
          double temp = enorm(neqs, wa3);
          double prered = 0.0;
          if( temp < fnorm ){
            prered = 1.0 - pow(temp/fnorm,2.0);
          }

          // compute the ratio of the actual to the predicted reduction.

          double ratio = 0.0;
          if( prered > 0.0 ){
            ratio = actred / prered;
          }

          // update the step bound.

          if( ratio < 0.1 ){
            ncsuc = 0;
            ncfail += 1;
            delta = 0.5 * delta;
          } else {
            ncfail = 0;
            ncsuc += 1;
            if( ratio >= 0.5 || ncsuc > 1 ){
              delta = fmax(delta, pnorm / 0.5);
            }
            if( fabs(ratio-1.0) <= 0.1 ){
              delta = pnorm / 0.5;
            }
          }

          // test for successful iteration.

          if( ratio >= 1.e-4 ){

            // successful iteration. update x, fvec, and their norms.

            for( int j=0 ; j<neqs ; ++j ){
              x[j] = wa2[j];
              wa2[j] = diag[j] * x[j];
              f[j] = wa4[j];
            }

            xnorm = enorm(neqs, wa2);
            fnorm = fnorm1;
            iter++;
          }

          // determine the progress of the iteration.

          nslow1++;
          if( actred >= 1.e-3 ){
            nslow1 = 0;
          }
          if( jeval ){
            nslow2++;
          }
          if( actred >= 0.1 ){
            nslow2 = 0;
          }

          // test for convergence.

          if( delta <= xtol * xnorm || fnorm == 0.0 ){
            info = 1;
          }
          if( info != 0 ){
            finished = 1;
            break;
          }

          // tests for termination and stringent tolerances.

          if( nfev >= maxfev ){
            info = 2;
          }
          if( 0.1 * fmax(0.1*delta, pnorm) <= epsmch * xnorm ){
            info = 3;
          }
          if( nslow2 == 5 ){
            info = 4;
          }
          if( nslow1 == 10 ){
            info = 5;
          }

          if( info != 0 ){
            finished = 1;
            break;
          }

          // criterion for recalculating jacobian.

          if( ncfail == 2 ){
            break;
          }

          // calculate the rank one modification to the jacobian
          // and update qtf if necessary.
          if( pnorm == 0.0 ){
            break;
          }

          for( int j=0 ; j<neqs ; ++j ){
            double sum = 0.0;
            for( int i=0 ; i<neqs ; ++i ){
              sum += jac[i][j] * wa4[i];
            }
            wa2[j] = (sum - wa3[j]) / pnorm;
            wa1[j] = diag[j] * ((diag[j] * wa1[j]) / pnorm);
            if( ratio >= 1.e-4 ){
              qtf[j] = sum;
            }
          }

          // compute the qr factorization of the updated jacobian.

          r1updt(neqs, r, wa1, wa2, wa3, &sing);
          r1mpyq_nxn(neqs, jac, wa2, wa3);
          r1mpyq(neqs, qtf, wa2, wa3);

          // end of the inner loop.

          jeval = 0;
        }

        // end of the outer loop

        if( finished == 1 ){
          break;
        }
      }

      // termination, either normal or user imposed.

      if( iflag < 0 ){
        info = iflag;
      }
      iflag = 0;

////////////////////////////////////////////////////////////////////////////////////
      
      fcn(x, f, state);

      // Check for convergence

      int converged = 1; //true

      if( neqs == 3 ){
        if( x[2] < 1.0 || x[2] >= 20.0 ){
          // If temperature is too low or even negative.
          // Or if the temperature is too high
          // Then the solved Temperature is not correct.
          converged = 0;
        }
      }

      for( int n=0 ; n<neqs ; ++n ){
        if( fabs(f[n]) >= eps) {
          converged = 0;
          break;
        }
      }

      if( converged == 1 ){
        state->mu_p = x[0];
        state->mu_n = x[1];
        if( neqs == 3 ){
          // Convert to T from T9
          state->T = x[2] * 1.e9;
        }
        return;
      }

      is_pos_new = f[0] > 0.0 && f[1] > 0.0;

      if( is_pos_old != is_pos_new ){ dx *= 0.8; }

      if( is_pos_new == 1 ){
        inner_x[0] -= dx;
        inner_x[1] -= dx;
      }else{
        inner_x[0] += dx;
        inner_x[1] += dx;
      }

      is_pos_old = is_pos_new;

    }

    outer_x[0] -= 0.5;

  }

  // if (hj.info != 1) {
  //     amrex::Error("failed to solve");
  // }
//#ifndef AMREX_USE_GPU
//    std::cout << "NSE solver failed with these conditions: " << std::endl;
//    std::cout << "Temperature: " << state_data.state.T << std::endl;
//    std::cout << "Density: " << state_data.state.rho << std::endl;
//    std::cout << "Ye: " << state_data.state.y_e << std::endl;
//    std::cout << "Internal Energy: " << state_data.state.e << std::endl;
//    std::cout << "Initial mu_p and mu_n: " << state_data.state.mu_p
//              << ", " << state_data.state.mu_n << std::endl;
//#endif

  printf("failed to solve\n");
  exit(0);
}

// A newton-raphson solver for finding nse state used for calibrating
// chemical potential of proton and neutron
void nse_nr_solver( int neqs , struct NSE * state , double eps ){
  // Currently NR solver only works for (rho, T, Ye) input

  if( neqs != 2 ){
    printf("NR solver currently doesn't solve for (rho, e, Ye) input. Please use hybrid solver with nse.use_hybrid_solver=1.");
    exit(0);
  }

  // whether nse solver converged or not
  int converged = 0;

  double f[2];
  double jac[2][2];
  double x[2];

  x[0] = state->mu_p;
  x[1] = state->mu_n;

  jcn(x, jac, state);
  fcn(x, f, state);

  // store determinant for finding inverse jac
  double det;
  // store inverse jacobian
  double inverse_jac[2][2];
  // difference in chemical potential of proton and neutron
  double d_mu_p = 1e10;
  double d_mu_n = 1e10;

  // begin newton-raphson
  for( int i=0 ; i<500 ; ++i ){

    // check if current state fulfills constraint equation
    if( fabs(d_mu_p) < eps * fabs(x[0]) && fabs(d_mu_n) < eps * fabs(x[1]) ){
      converged = 1;
      state->mu_p = x[0];
      state->mu_n = x[1];
      printf("%d\n",i);
      break;
    }

    // Find the max of the jacobian used for scaling determinant to prevent digit overflow
    double scale_fac = fmax(jac[1][1],fmax(jac[1][0],fmax(jac[0][0],jac[0][1])));

    // if jacobians are small, then no need for scaling
    if( scale_fac < 1.0e150 ){ scale_fac = 1.0; }

    // Specific inverse 2x2 matrix, perhaps can write a function for solving n systems of equations.
    det = jac[0][0] / scale_fac * jac[1][1] - jac[0][1] / scale_fac * jac[1][0];

    // check if determinant is 0
    if( det == 0.0 ){ 
      printf("Jacobian is a singular matrix! Try a different initial guess!"); 
      exit(0);
    }

    // find inverse jacobian
    inverse_jac[0][0] = jac[1][1] / scale_fac / det;
    inverse_jac[0][1] = -jac[0][1] / scale_fac / det;
    inverse_jac[1][0] = -jac[1][0] / scale_fac / det;
    inverse_jac[1][1] = jac[0][0] / scale_fac / det;

    // find the difference
    d_mu_p = -(f[0] * inverse_jac[0][0] + f[1] * inverse_jac[0][1]);
    d_mu_n = -(f[0] * inverse_jac[1][0] + f[1] * inverse_jac[1][1]);

    // if diff goes beyond 1.0e3, likely that its not making good progress..
    if( fabs(d_mu_p) > 1.0e3 || fabs(d_mu_n) > 1.0e3 ){ printf("Not making good progress, breaking"); }

    // update new solution
    x[0] += d_mu_p;
    x[1] += d_mu_n;

    // check whether solution results in nan
    if( isnan(x[0]) || isnan(x[1]) ){ printf("Nan encountered, likely due to overflow in digits or not making good progress"); }

    // update constraint
    jcn(x, jac, state);
    fcn(x, f, state);
  }

  if( !converged ){ printf("NSE solver failed to converge!"); }

}

void nse_solve( int neqs , struct NSE * nse_state , double eps ){
  // Interface for calling the Hybrid-Powell or Netwon-Raphson Solver
  int singular_network;
  int use_hybrid_solver = 1;

  if( use_hybrid_solver == 1 ){
    nse_hybrid_solver( neqs , nse_state , eps );
  } else {
    singular_network = 1;
    for( int n=0 ; n<NUM_I ; ++n ){
      if( zion[n] != aion[n] - zion[n] ){
        singular_network = 0;
      }
    }

    if (singular_network) {
      printf("This network always results in singular jacobian matrix, thus can't find nse mass fraction using nr!");
      exit(0);
    }

    nse_nr_solver( neqs , nse_state , eps );
  }
}

///////////////////////////////////////////////////////////////////////////////////

void get_nse_state_rty( struct NSE * state, double eps ){
  // Find the NSE state using (rho, T, Ye) as input.
  // This returns NSE mass fraction and
  // solved chemical potentials.
  // It also updates chemical potentials to the input state.
 
  // First get the temperature-independent part.
  get_T_independent_nse_state( state );

  // Since Temperature is fixed in this mode,
  // Get Temperature dependent exponent pieces, including screening
  // The values are stored in log_data
  get_T_dependent_nse_state( state );

  // invoke newton-raphson or hybrj to solve chemical potential of proton and neutron
  // which are the exponent part of the nse mass fraction calculation
  int neqs = 2;
  nse_solve( neqs , state , eps );

  // Apply the log_data for calculating nse mass fractions
  // This updates y_e as well.
  apply_nse_exponent( state );

  // update rhoX in the output nse_state
  for( int n=0 ; n<NUM_I ; ++n ){
    state->yn[n] = state->rho * state->xn_exp[n];
  }
}

double nse_abar( double T, double rho, double Ye, double mu_p, double mu_n, double * comp ){
  ///
  /// This function calculates abar from NSE using
  /// Temp, rho, and Ye
  /// mu_p and mu_n are initial guesses

  struct NSE nse_state;
  nse_state.rho = rho;
  nse_state.y_e = Ye;
  nse_state.y_e_init = Ye;
  nse_state.T = T;
  nse_state.mu_p = mu_p;
  nse_state.mu_n = mu_n;

  get_nse_state_rty( &nse_state, 1.0e-10 );

  double abar = 0.0;
  for( int n=0 ; n<NUM_I ; ++n ){
    abar += nse_state.xn[n] * aion_inv[n];
    comp[n] = nse_state.xn[n];
  }
  abar = 1.0 / abar;

  return abar;
}

double nse_dabar_dT( double T, double rho, double Ye, double mu_p, double mu_n, double * dXdT ){
  ///
  /// This function constructs dabar_dT
  /// This should be 2nd order accurate
  ///

  // deviation in temperature
  double dT = 1.0e-6 * T;

  // Calculate derivative using five-point stencil method
  // double dabar_dT = (-nse_abar(T + 2.0_rt*dT, rho, Ye, mu_p, mu_n) +
  //                          8.0_rt * nse_abar(T+dT, rho, Ye, mu_p, mu_n) -
  //                          8.0_rt * nse_abar(T-dT, rho, Ye, mu_p, mu_n) +
  //                          nse_abar(T - 2.0_rt*dT, rho, Ye, mu_p, mu_n)
  //                         ) / (12.0_rt * dT);

  // Calculate derivative using central differencing
  double X_p[NUM_I], X_m[NUM_I];
  double a_p = nse_abar(T + dT, rho, Ye, mu_p, mu_n, X_p);
  double a_m = nse_abar(T - dT, rho, Ye, mu_p, mu_n, X_m);

  double dabar_dT = 0.5 * ( a_p - a_m ) / dT;
  //for( int n=0 ; n<NUM_I ; ++n ){ dXdT[n] = 0.5 * ( X_p[n] - X_m[n] ) / dT; }
  double ebind_p = 0.0, ebind_m = 0.0;
  for( int n=0 ; n<NUM_I ; ++n ){ 
    ebind_p += EBIND[n]*X_p[n];
    ebind_m += EBIND[n]*X_m[n];
  }
  dXdT[0] = 0.5 * ( ebind_p - ebind_m ) / dT;

  return dabar_dT;
}

void nse_T_from_e( double rho, double e_in, double Ye, double * T, double * mu_p, double * mu_n ){
  ///
  /// This function inverts this form of the EOS to find the T
  /// that satisfies the EOS and NSE given an input e and rho.
  ///
  /// if we are in NSE, then the entire thermodynamic state is just
  /// a function of rho, T, Ye.  We can write the energy as:
  ///
  ///    e = e(rho, T, Y_e, Abar(rho, T, Ye))
  ///
  /// where we note that Abar is a function of those same inputs.
  ///
  /// The basic idea is that Abar and Zbar are both functions of
  /// rho, T, Ye through NSE calculations, so we express the energy
  /// as:
  ///
  ///      e = e(rho, T, Abar(rho, T, Ye), Zbar(rho, T, Ye)
  ///
  /// and NR on that.  Note that Zbar = Ye Abar, so we can group
  /// those derivative terms together.
  ///
  /// T come in as initial guesses and are updated on output
  ///

  double ttol = 1e-8;
  int max_iter = 100;

  int converged = 0; //false

  int iter = 0;

  double temp = *T;
  // initialize burn_state
  struct NSE nse_state;
  nse_state.rho = rho;
  nse_state.y_e = Ye;
  nse_state.y_e_init = Ye;
  nse_state.mu_p = *mu_p;
  nse_state.mu_n = *mu_n;

  double f,df,dT,e_eos,e_bind;
  double prim[NUM_Q],dXdT[NUM_I];
  double derivs[3],dedT,dedA,dedZ;

  while( converged == 0 && iter < max_iter ){

    // update Temperature
    nse_state.T = temp;

    get_nse_state_rty( &nse_state , 1.0e-10 );

    // Call the EOS to get internal energy
    prim[RHO] = rho;
    for( int n=0 ; n<NUM_I ; ++n ){ prim[XXX+n] = nse_state.xn_exp[n]; }
    e_eos = get_eint( prim , &temp );
    get_Ederivs( prim , &temp , derivs );
    dedT = derivs[0];
    dedA = derivs[1];
    dedZ = derivs[2];

    double dabar_dT = nse_dabar_dT(temp, nse_state.rho, nse_state.y_e, nse_state.mu_p, nse_state.mu_n, dXdT);
    
    // f is the quantity we want to zero
    f = e_eos - e_in;
    df = dedT + dedA * dabar_dT + nse_state.y_e * dedZ * dabar_dT;
    //printf("df_eos = %e\n",df);
    if( solve_nse_e_mode == 2 ){
      e_bind = 0.0;
      for( int n=0 ; n<NUM_I ; ++n ){
        e_bind += EBIND[n]*nse_state.xn_exp[n];
        //df += EBIND[n]*dXdT[n];
        //printf("dXdT = %e\n",dXdT[n]);
      }
      f += e_bind;
      df += dXdT[0];
    }
    //printf("df_tot = %e\n",df);
    //printf("eint: %e , ebind: %e -> etot: %e\n", e_eos, e_bind, e_eos+e_bind);
    // compute the correction to our guess
    dT = -f / df;
    //printf("dT = %e\n",dT);

    // update the temperature
    temp = fmin(fmax(temp + dT, 0.25*temp), 4.0*temp);

    // check convergence
    if( fabs(dT) < ttol*temp ){
      converged = 1; //true
      *T = temp;
      *mu_p = nse_state.mu_p;
      *mu_n = nse_state.mu_n;
      printf("final Temp = %e\n",temp);
    }
    iter++;
  }
}

void get_nse_state_rey( struct NSE * state, double eps ){
  // Find the NSE state using (rho, e, Ye) as input.
  // We can solve in two ways:
  // 1. Find the NSE state using (rho, T, Ye), then do NR iterations
  // on that to find the correct T that matches e.
  // 2. Find the NSE state directly using (rho, e, Ye),
  // where now we have an additional constraint equation on e.
  // So there are three constraint equations and a 3x3 Jacobian.

  if( solve_nse_e_mode == 1 ){
    // Option 1

    // Get T0 such that it is consistent with e.
    // This also updates mu_p and mu_n of the input state.

    double T0 = state->T;
    nse_T_from_e(state->rho, state->e, state->y_e, &T0, &state->mu_p, &state->mu_n);

    // Create burn state with the updated Temperature
    // Then find NSE state via (rho, T, Ye)

    state->T = T0;
    get_nse_state_rty( state , eps );

    return;
  } else if( solve_nse_e_mode == 2 ){
    // Option 2

    // Get T0 such that it is consistent with etot = eint+ebind.
    // This also updates mu_p and mu_n of the input state.

    double T0 = state->T;
    nse_T_from_e(state->rho, state->etot, state->y_e, &T0, &state->mu_p, &state->mu_n);

    // Create burn state with the updated Temperature
    // Then find NSE state via (rho, T, Ye)

    state->T = T0;
    get_nse_state_rty( state , eps );

    return;
  }

  // Option 3

  // First get the temperature-independent part.
  get_T_independent_nse_state( state );

  // invoke newton-raphson or hybrj to solve mu_p, mu_n and T.
  // For (rho, e, Ye) input, there are three constraint eqns.
  // Note that state_data.state already pre-computed T-independent terms.

  int neqs = 3;
  nse_solve( neqs , state , eps );

  // Now that mu_p, mu_n and T are solved and stored in state_data.state.
  // Finish evaluating the NSE composition.

  // Evaluate the temperature dependent piece, including screening
  get_T_dependent_nse_state( state );

  // Apply exponent for calculating NSE mass fractions
  apply_nse_exponent( state );

  // update rhoX in the output nse_state
  for( int n=0 ; n<NUM_I ; ++n ){
    state->yn[n] = state->rho * state->xn_exp[n];
  }

  return;
}

// Get the NSE state;
void get_actual_nse_state( double y[] , double rho , double T , double eps ){

  // Convert molar abundance to mass density for the NSE solver
  for( int n=0 ; n<NUM_I ; ++n ){ y[n] = rho * y[n] * aion[n]; }

  // Find input electron fraction
  double y_e = 0.0;
  for( int n=0 ; n<NUM_I ; ++n ){ y_e += y[n] * zion[n] * aion_inv[n]; }
  y_e /= rho;

  struct NSE state;
  state.rho = rho;
  state.T = T;
  for( int n=0 ; n<NUM_I ; ++n ){ state.yn[n] = y[n]; }
  state.y_e = y_e;
  state.y_e_init = y_e;
  state.mu_p = -5.;//-7.0;  //initial guess
  state.mu_n = -10.;//-11.0;

  double prim[NUM_Q],e_bind;

  if (input == 0) {
    // Using (rho, T, Ye) as input
    get_nse_state_rty( &state , eps );
  } else if (input == 1) {
    // Using (rho, e, Ye) as input
    prim[RHO] = state.rho;
    e_bind = 0.0;
    for( int n=0 ; n<NUM_I ; ++n ){ 
      prim[XXX+n] = state.yn[n]/state.rho;
      e_bind += EBIND[n]*prim[XXX+n];
    }
    state.e = get_eint( prim , &state.T );
    state.etot = state.e + e_bind;
    printf("eint: %e , ebind: %e -> etot: %e\n", state.e, e_bind, state.etot);
    get_nse_state_rey( &state , eps );
  } else {
    // When getting an unknown input
    printf("Unknown NSE input mode");
    exit(0);
  }

  // Output //
  for( int n=0 ; n<NUM_I ; ++n ){ y[n] = state.yn[n] * aion_inv[n] / state.rho; }

}
