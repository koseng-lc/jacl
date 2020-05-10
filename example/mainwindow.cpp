/**
*   @author : koseng (Lintang)
*   @brief : GUI Source
*/

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    //-- DC motor parameters - Maxon 2322.983-11.225-200
    , bm(0.0436e-3), jm(5.71e-7), ki(30.7e-3), la(1.97e-3), kb(30.7e-3), ra(21.6)
    , weight_({0.1, 0.1, 0.1, 0.1, 0.1, 0.1})
    //-- MIMO DC motor
    , ss_(&bm, &jm, &ki, &la, &kb, &ra)    
    , csys_(&ss_)
    , csys_plt_(&csys_,{1,2,3,4,5,6,7,8})
    , cobserver_(&ss_, arma::zeros<arma::mat>(3,3))
    , cobserver_plt_(&cobserver_, {1,2,3,6,7,8})
    //-- SIMO DC motor
    , simo_(&bm, &jm, &ki, &la, &kb, &ra)
    , dsys_simo_(&dsimo_)
    , dobserver_simo_(&dsimo_, arma::zeros<arma::mat>(3,3))
    , dsys_simo_plt_(&dsys_simo_, {1,2,3,4}, SAMPLING_PERIOD)
    , dobserver_simo_plt_(&dobserver_simo_, {1,2,3}, SAMPLING_PERIOD)
    , ifd_(&dsys_simo_, {10.,9.,8.})
    , sifd_(&dsys_simo_, {10., .05, .07})
    , m_real_(&bm,&jm,&ki,&la,&kb,&ra)
    , m_sys_(&m_icm_, SAMPLING_PERIOD)
    , pos_sys_(&pos_icm_, SAMPLING_PERIOD)
    , pos_dctrl_sys_(&pos_dctrl_, SAMPLING_PERIOD)
    , pos_dctrl_plt_(&pos_dctrl_sys_, {4,5}, SAMPLING_PERIOD)
    //--
    , G_(&bm, &jm, &ki, &la, &kb, &ra)
    , ref_(3, 1, arma::fill::zeros)
    , cl_status_(true)
    , posctrl_sys_(&k_pos_)
    , posctrl_plt_(&posctrl_sys_,{7,8})
    , spdctrl_sys_(&k_spd_)
    , spdctrl_plt_(&spdctrl_sys_,{6,7})
    , pg_(9.8), pl_(1.0), pk_(0.1), pm_(0.5)
    , nlp_(&pg_, &pl_, &pk_, &pm_)
    , nlp_sys_(&nlp_){
    //-- GUI    
    ui->setupUi(this);
    this->removeToolBar(ui->mainToolBar);
    QFile style_file("../example/gui/dark_style_lintang.qss");
    style_file.open(QFile::ReadOnly);
    QString style_sheet(style_file.readAll());
    this->setStyleSheet(style_sheet);

    //-- DC Motor Open-Loop
    std::cout << "Preparing system ..." << std::endl;
    {
        LinearStateSpace::Formula A11 = JC(ss_,.0); LinearStateSpace::Formula A12 = JC(ss_,1.0); LinearStateSpace::Formula A13 = JC(ss_,.0);
        LinearStateSpace::Formula A21 = JC(ss_,.0); LinearStateSpace::Formula A22 = JE(ss_,-_p(iBm)/_p(iJm)); LinearStateSpace::Formula A23 = JE(ss_,_p(iKi)/_p(iJm));
        LinearStateSpace::Formula A31 = JC(ss_,.0); LinearStateSpace::Formula A32 = JE(ss_,-_p(iKb)/_p(iLa)); LinearStateSpace::Formula A33 = JE(ss_,-_p(iRa)/_p(iLa));

        LinearStateSpace::Formulas fA{
            A11,A12,A13,
            A21,A22,A23,
            A31,A32,A33
        };

        LinearStateSpace::Formula B11 = JC(ss_, .0); LinearStateSpace::Formula B12 = JC(ss_, .0);
        LinearStateSpace::Formula B21 = JC(ss_, .0); LinearStateSpace::Formula B22 = JE(ss_, -1.0/_p(iJm));
        LinearStateSpace::Formula B31 = JE(ss_, 1.0/_p(iLa)); LinearStateSpace::Formula B32 = JC(ss_, .0);

        LinearStateSpace::Formulas fB{
            B11, B12,
            B21, B22,
            B31, B32
        };

        LinearStateSpace::Formula C11 = JC(ss_, 1.0); LinearStateSpace::Formula C12 = JC(ss_, .0); LinearStateSpace::Formula C13 = JC(ss_, .0);
        LinearStateSpace::Formula C21 = JC(ss_, .0); LinearStateSpace::Formula C22 = JC(ss_, 1.0); LinearStateSpace::Formula C23 = JC(ss_, .0);
        LinearStateSpace::Formula C31 = JC(ss_, .0); LinearStateSpace::Formula C32 = JC(ss_, .0); LinearStateSpace::Formula C33 = JC(ss_, 1.0);

        LinearStateSpace::Formulas fC{
            C11, C12, C13,
            C21, C22, C23,
            C31, C32, C33
        };

        LinearStateSpace::Formula D11 = JC(ss_, .0); LinearStateSpace::Formula D12 = JC(ss_, .0);
        LinearStateSpace::Formula D21 = JC(ss_, .0); LinearStateSpace::Formula D22 = JC(ss_, .0);
        LinearStateSpace::Formula D31 = JC(ss_, .0); LinearStateSpace::Formula D32 = JC(ss_, .0);

        LinearStateSpace::Formulas fD{
            D11, D12,
            D21, D22,
            D31, D32,
        };

        ss_.setA(fA);
        ss_.setB(fB);
        ss_.setC(fC);
        ss_.setD(fD);
    }

    ss_.A().print("A : ");
    ss_.B().print("B : ");
    ss_.C().print("C : ");
    ss_.D().print("D : ");

    std::cout << "Controllable : " << jacl::common::controllable(ss_.A(), ss_.B()) << std::endl;
    std::cout << "Observable : " << jacl::common::observable(ss_.A(), ss_.C()) << std::endl;

    //-- Observer for MIMO DC motor
    arma::mat observer_K;
    arma::mat cpoles{-50,-51,-52};
    jacl::pole_placement::KautskyNichols(&ss_, cpoles, &observer_K, jacl::pole_placement::PolePlacementType::Observer);
    cobserver_.setGain(observer_K.t());

    //-- G Realization
    {
        GRealization::Formula A11, A12, A13;
        A11 = JC(G_, .0); A12 = JC(G_, 1.0); A13 = JC(G_, .0);
        GRealization::Formula A21, A22, A23;
        A21 = JC(G_, .0); A22 = JE(G_, -_p(iBm)/_p(iJm)); A23 = JE(G_, _p(iKi)/_p(iJm));
        GRealization::Formula A31, A32, A33;
        A31 = JC(G_, .0); A32 = JE(G_, -_p(iKb)/_p(iLa)); A33 = JE(G_, -_p(iRa)/_p(iLa));

        GRealization::Formula B11, B12, B13, B14, B15, B16, B17, B18;
        B11 = B12 = B13 = B14 = B15 = B16 = B17 = B18 = JC(G_, .0);
        GRealization::Formula B21, B22, B23, B24, B25, B26, B27, B28;
        B21 = B22 = B23 = JE(G_, -1.0/_p(iJm)); B24 = B25 = B26 = B27 = JC(G_, .0); B28 = JE(G_, -1.0/_p(iJm));
        GRealization::Formula B31, B32, B33, B34, B35, B36, B37, B38;
        B31 = B32 = B33 = JC(G_, .0); B34 = B35 = B36 = B37 = JE(G_, -1.0/_p(iLa)); B38 = JC(G_, .0);

        GRealization::Formula C11, C12, C13;
        C11 = JC(G_, .0); C12 = JC(G_, 1.0); C13 = JC(G_, .0);
        GRealization::Formula C21, C22, C23;
        C21 = JC(G_, .0); C22 = JE(G_, -_p(iBm)/_p(iJm)); C23 = JE(G_, _p(iKi)/_p(iJm));
        GRealization::Formula C31, C32, C33;
        C31 = JC(G_, .0); C32 = JC(G_, .0); C33 = JC(G_, 1.0);
        GRealization::Formula C41, C42, C43;
        C41 = JC(G_, .0); C42 = JE(G_, -_p(iKb)/_p(iLa)); C43 = JE(G_, -_p(iRa)/_p(iLa));
        GRealization::Formula C51, C52, C53;
        C51 = JC(G_, .0); C52 = JC(G_, 1.0); C53 = JC(G_, .0);
        GRealization::Formula C61, C62, C63;
        C61 = JC(G_, .0); C62 = JC(G_, .0); C63 = JC(G_, 1.0);
        GRealization::Formula C71, C72, C73;
        C71 = JC(G_, 1.0); C72 = JC(G_, .0); C73 = JC(G_, .0);
        GRealization::Formula C81, C82, C83;
        C81 = JC(G_, .0); C82 = JC(G_, 1.0); C83 = JC(G_, .0);
        GRealization::Formula C91, C92, C93;
        C91 = JC(G_, .0); C92 = JC(G_, .0); C93 = JC(G_, 1.0);

        GRealization::Formula D11, D12, D13, D14, D15, D16, D17, D18;
        D11 = D12 = D13 = D14 = D15 = D16 = D17 = D18 = JC(G_, .0);
        GRealization::Formula D21, D22, D23, D24, D25, D26, D27, D28;
        D21 = D22 = D23 = JE(G_, -1/_p(iJm)); D24 = D25 = D26 = JC(G_, .0); D27 = JE(G_, -1/_p(iJm)); D28 = JC(G_, .0);
        GRealization::Formula D31, D32, D33, D34, D35, D36, D37, D38;
        D31 = D32 = D33 = D34 = D35 = D36 = D37 = D38 = JC(G_, .0);
        GRealization::Formula D41, D42, D43, D44, D45, D46, D47, D48;
        D41 = D42 = D43 = JC(G_, .0); D44 = D45 = D46 = JE(G_, -1/_p(iLa)); D47 = JC(G_, .0); D48 = JE(G_, -1/_p(iLa));
        GRealization::Formula D51, D52, D53, D54, D55, D56, D57, D58;
        D51 = D52 = D53 = D54 = D55 = D56 = D57 = D58 = JC(G_, .0);
        GRealization::Formula D61, D62, D63, D64, D65, D66, D67, D68;
        D61 = D62 = D63 = D64 = D65 = D66 = D67 = D68 = JC(G_, .0);
        GRealization::Formula D71, D72, D73, D74, D75, D76, D77, D78;
        D71 = D72 = D73 = D74 = D75 = D76 = D77 = D78 = JC(G_, .0);
        GRealization::Formula D81, D82, D83, D84, D85, D86, D87, D88;
        D81 = D82 = D83 = D84 = D85 = D86 = D87 = D88 = JC(G_, .0);
        GRealization::Formula D91, D92, D93, D94, D95, D96, D97, D98;
        D91 = D92 = D93 = D94 = D95 = D96 = D97 = D98 = JC(G_, .0);

        GRealization::Formulas fA{
            A11, A12, A13,
            A21, A22, A23,
            A31, A32, A33
        };

        GRealization::Formulas fB{
            B11, B12, B13, B14, B15, B16, B17, B18,
            B21, B22, B23, B24, B25, B26, B27, B28,
            B31, B32, B33, B34, B35, B36, B37, B38
        };

        GRealization::Formulas fC{
            C11, C12, C13,
            C21, C22, C23,
            C31, C32, C33,
            C41, C42, C43,
            C51, C52, C53,
            C61, C62, C63,
            C71, C72, C73,
            C81, C82, C83,
            C91, C92, C93,
        };

        GRealization::Formulas fD{
            D11, D12, D13, D14, D15, D16, D17, D18,
            D21, D22, D23, D24, D25, D26, D27, D28,
            D31, D32, D33, D34, D35, D36, D37, D38,
            D41, D42, D43, D44, D45, D46, D47, D48,
            D51, D52, D53, D54, D55, D56, D57, D58,
            D61, D62, D63, D64, D65, D66, D67, D68,
            D71, D72, D73, D74, D75, D76, D77, D78,
            D81, D82, D83, D84, D85, D86, D87, D88,
            D91, D92, D93, D94, D95, D96, D97, D98,
        };

        G_.setA(fA); // use A from Plant
        G_.setB(fB);
        G_.setC(fC);
        G_.setD(fD);
    }

//    G_.A().print("\nGA : ");
//    G_.B().print("GB : ");
//    G_.C().print("GC : ");
//    G_.D().print("GD : ");

    //-- P-Delta Realization
    {
        PRealization::Formulas fA{
            JC(P_, -1790), JC(P_, 945.2), JC(P_, -0.4491), JC(P_, 9.531), JC(P_, -0.6474), JC(P_, 0.02418), JC(P_, -4.527e-05), JC(P_, 4.543e-08), JC(P_, 8.589e-15),
            JC(P_, 2266), JC(P_, -1537), JC(P_, -0.5102), JC(P_, 10.96), JC(P_, -0.7445), JC(P_, 0.02781), JC(P_, -5.206e-05), JC(P_, 5.225e-08), JC(P_, -1.226e-14),
            JC(P_, -0.238), JC(P_, 0.313), JC(P_, -3136), JC(P_, 80.45), JC(P_, -5.465), JC(P_, 0.2041), JC(P_, -0.0003821), JC(P_, 3.835e-07), JC(P_, -3.118e-14),
            JC(P_, -1.281), JC(P_, -1.065), JC(P_, 0.4703), JC(P_, -3128), JC(P_, 207.8), JC(P_, 0.009844), JC(P_, -1.843e-05), JC(P_, 1.85e-08), JC(P_, -5.121e-14),
            JC(P_, 0.08699), JC(P_, 0.07231), JC(P_, -0.03195), JC(P_, 207.8), JC(P_, -77.96), JC(P_, 110.8), JC(P_, 1.252e-06), JC(P_, -1.256e-09), JC(P_, -1.43e-14),
            JC(P_, -0.003249), JC(P_, -0.002701), JC(P_, 0.001193), JC(P_, -0.4271), JC(P_, 104.3), JC(P_, -180.3), JC(P_, 5.586), JC(P_, 4.693e-11), JC(P_, -1.897e-14),
            JC(P_, 6.083e-06), JC(P_, 5.057e-06), JC(P_, -2.234e-06), JC(P_, 13.67), JC(P_, 189), JC(P_, -322.2), JC(P_, -161.8), JC(P_, 2.968), JC(P_, -7.243e-14),
            JC(P_, -6.105e-09), JC(P_, -5.075e-09), JC(P_, 2.242e-09), JC(P_, -6.387), JC(P_, -88.26), JC(P_, 154.4), JC(P_, -12.19), JC(P_, -188.8), JC(P_, 3.215e-15),
            JC(P_, 1.277e-14), JC(P_, -7.597e-15), JC(P_, -5.757e-16), JC(P_, -0.05278), JC(P_, -0.7943), JC(P_, -0.4598), JC(P_, -0.01586), JC(P_, -0.0002494), JC(P_, -187.4)
        };

        PRealization::Formulas fB{
            JC(P_, 0.6773), JC(P_, -4.319e-18),
            JC(P_, -0.4418), JC(P_, 2.803e-18),
            JC(P_, 3.676), JC(P_, -1.545e-17),
            JC(P_, -63.75), JC(P_, -0.001311),
            JC(P_, 4.33), JC(P_, -0.01824),
            JC(P_, -0.1617), JC(P_, 0.02794),
            JC(P_, 0.0003028), JC(P_, -0.1938),
            JC(P_, -3.039e-07), JC(P_, -0.4921),
            JC(P_, -1.017e-15), JC(P_, -1.313)
        };

        PRealization::Formulas fC{
            JC(P_, 0.0371), JC(P_, 0.03426), JC(P_, -7.107e-06), JC(P_, 0.05756), JC(P_, 0.8638), JC(P_, 0.5002), JC(P_, 0.01814), JC(P_, 0.002169), JC(P_, -0.004907),
            JC(P_, -2.435e-05), JC(P_, 1.049e-05), JC(P_, 0.6569), JC(P_, 0.03607), JC(P_, -0.0252), JC(P_, 0.02546), JC(P_, 0.3908), JC(P_, 0.8316), JC(P_, 0.3932),
            JC(P_, -0.5196), JC(P_, 0.3392), JC(P_, -2.416), JC(P_, 50.2), JC(P_, -3.427), JC(P_, 0.1587), JC(P_, -0.05068), JC(P_, 0.02228), JC(P_, -1.364e-15)
        };

        PRealization::Formulas fD{
            JC(P_, .0), JC(P_, .0),
            JC(P_, .0), JC(P_, .0),
            JC(P_, .0), JC(P_, .0)
        };

        P_.setA(fA);
        P_.setB(fB);
        P_.setC(fC);
        P_.setD(fD);
    }

    //-- Realization of Inter-connection Matrix    
    {
        arma::mat temp1, temp2, temp3;
        arma::mat zeros9x3(9,3,arma::fill::zeros);
        arma::mat zeros2x9(2,9,arma::fill::zeros);
        arma::mat zeros3x3(3,3,arma::fill::zeros);
        arma::mat zeros3x2(3,2,arma::fill::zeros);
        arma::mat zeros2x3(2,3,arma::fill::zeros);
        arma::mat zeros2x2(2,2,arma::fill::zeros);
        arma::mat eye3x3(3,3,arma::fill::eye);
        arma::mat eye2x2(2,2,arma::fill::eye);

        temp1 = arma::join_horiz(zeros9x3, P_.B());
        temp2 = arma::join_horiz(zeros9x3, temp1);
        arma::mat B = arma::join_horiz(P_.B(), temp2);

        temp1 = arma::join_vert(zeros2x9, -1.*P_.C());
        arma::mat C = arma::join_vert(-1.*P_.C(), temp1);

        temp1 = arma::join_horiz(zeros3x3, eye3x3);
        temp2 = arma::join_horiz(zeros3x2, temp1);
        temp3 = arma::join_horiz(zeros2x2, arma::join_horiz(zeros2x3, zeros2x3));
        arma::mat D11 = arma::join_vert(temp2, temp3);

        arma::mat D12 = arma::join_vert(zeros3x2, eye2x2);

        temp1 = arma::join_horiz(-1.*eye3x3,eye3x3);
        arma::mat D21 = arma::join_horiz(zeros3x2,temp1);
        
        arma::mat D22 = zeros3x2;

        temp1 = arma::join_horiz(D11, D12);
        temp2 = arma::join_horiz(D21, D22);
        arma::mat D = arma::join_vert(temp1, temp2);

        ICM_.setA(P_.A());
        ICM_.setB(B);
        ICM_.setC(C);
        ICM_.setD(D);

    }    

//    ICM_.A().print("ICM A : ");
//    ICM_.B().print("ICM B : ");
//    ICM_.C().print("ICM C : ");
//    ICM_.D().print("ICM D : ");

    // hinf_ = new HInf(&ICM_, 1.6204);
    // jacl::common::StateSpacePack K( hinf_->solve() );

    setupSIMODCMotor();
    setupPositionController();
    setupSpeedController();    

    //-- Simulator
    csys_plt_.init();
    csys_plt_.setTitle("DC Motor");
    csys_plt_.setDelay() = SAMPLING_PERIOD;
    csys_plt_.setPlotName({"Angular Position", "Angular Velocity", "Current",
                       "Voltage In", "Torque In",
                       "Angular Position", "Angular Velocity", "Current"});

    cobserver_plt_.init();
    cobserver_plt_.setTitle("Full-order Luenberger Observer of DC Motor");
    cobserver_plt_.setDelay() = SAMPLING_PERIOD;
    cobserver_plt_.setPlotName({"Est. Position", "Est. Velocity", "Est. Current"
                              ,"Est. Out Position", "Est. Out Velocity", "Est. Out Current"});

    setupWidgets();
    setupControllerDialog();
    setupMenus();
    setupActions();

    std::vector<double> y;
    auto dt(M_PI/180.);
    for(int i(0); i < 1000; i++){
        y.push_back(sin((double)i*dt));
    }
    jacl::plot(y,dt,"TEST",{"XY"});

    cl_thread_ = boost::thread{boost::bind(&MainWindow::closedLoopProcess, this)};

}

MainWindow::~MainWindow(){
    cl_status_ = false;
    cl_thread_.join();
    
    delete ui;
}

void MainWindow::setupWidgets(){

    //-- Params

    params_gb_ = new QGroupBox;

    params_gl_ = new QGridLayout;

    params_label_[iBm] = new QLabel;
    params_label_[iBm]->setText(tr("Viscous Friction : "));
    params_dsb_[iBm] = new QDoubleSpinBox;
    params_dsb_[iBm]->setSingleStep(1e-6);
    params_dsb_[iBm]->setDecimals(6);
    params_dsb_[iBm]->setMaximum(1.0);
    params_dsb_[iBm]->setMinimum(.0);
    params_dsb_[iBm]->setValue(bm.nominal);
    params_dsb_[iBm]->adjustSize();
//    params_dsb_[iBm]->resize(params_dsb_[iBm]->minimumSizeHint().width(),
//                params_dsb_[iBm]->minimumSizeHint().height());

    params_label_[iJm] = new QLabel;
    params_label_[iJm]->setText(tr("Rotor Inertia : "));
    params_dsb_[iJm] = new QDoubleSpinBox;
    params_dsb_[iJm]->setSingleStep(1e-7);
    params_dsb_[iJm]->setDecimals(7);
    params_dsb_[iJm]->setMaximum(1.0);
    params_dsb_[iJm]->setMinimum(.0);
    params_dsb_[iJm]->setValue(jm.nominal);
    params_dsb_[iJm]->adjustSize();

    params_label_[iKi] = new QLabel;
    params_label_[iKi]->setText(tr("Torque Constant : "));
    params_dsb_[iKi] = new QDoubleSpinBox;
    params_dsb_[iKi]->setSingleStep(1e-3);
    params_dsb_[iKi]->setDecimals(3);
    params_dsb_[iKi]->setMaximum(10.0);
    params_dsb_[iKi]->setMinimum(.0);
    params_dsb_[iKi]->setValue(ki.nominal);
    params_dsb_[iKi]->adjustSize();

    params_label_[iLa] = new QLabel;
    params_label_[iLa]->setText(tr("Armature Inductance : "));
    params_dsb_[iLa] = new QDoubleSpinBox;
    params_dsb_[iLa]->setSingleStep(1e-5);
    params_dsb_[iLa]->setDecimals(5);
    params_dsb_[iLa]->setMaximum(10.0);
    params_dsb_[iLa]->setMinimum(.0);
    params_dsb_[iLa]->setValue(la.nominal);
    params_dsb_[iLa]->adjustSize();

    params_label_[iKb] = new QLabel;
    params_label_[iKb]->setText(tr("Back EMF : "));
    params_dsb_[iKb] = new QDoubleSpinBox;
    params_dsb_[iKb]->setSingleStep(1e-3);
    params_dsb_[iKb]->setDecimals(3);
    params_dsb_[iKb]->setMaximum(10.0);
    params_dsb_[iKb]->setMinimum(.0);
    params_dsb_[iKb]->setValue(kb.nominal);
    params_dsb_[iKb]->adjustSize();

    params_label_[iRa] = new QLabel;
    params_label_[iRa]->setText(tr("Armature Resistance : "));
    params_dsb_[iRa] = new QDoubleSpinBox;
    params_dsb_[iRa]->setSingleStep(1e-2);
    params_dsb_[iRa]->setDecimals(2);
    params_dsb_[iRa]->setMaximum(50.0);
    params_dsb_[iRa]->setMinimum(.0);
    params_dsb_[iRa]->setValue(ra.nominal);
    params_dsb_[iRa]->adjustSize();

    params_gl_->addWidget(params_label_[iBm],0,0,1,1);
    params_gl_->addWidget(params_label_[iJm],1,0,1,1);
    params_gl_->addWidget(params_label_[iKi],2,0,1,1);
    params_gl_->addWidget(params_label_[iLa],3,0,1,1);
    params_gl_->addWidget(params_label_[iKb],4,0,1,1);
    params_gl_->addWidget(params_label_[iRa],5,0,1,1);
    params_gl_->addWidget(params_dsb_[iBm],  0,1,1,1);
    params_gl_->addWidget(params_dsb_[iJm],  1,1,1,1);
    params_gl_->addWidget(params_dsb_[iKi],  2,1,1,1);
    params_gl_->addWidget(params_dsb_[iLa],  3,1,1,1);
    params_gl_->addWidget(params_dsb_[iKb],  4,1,1,1);
    params_gl_->addWidget(params_dsb_[iRa],  5,1,1,1);

    params_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),6,2);

    params_gb_->setLayout(params_gl_);
    params_gb_->setTitle(tr("Parameters"));

    //-- Command

    command_gb_ = new QGroupBox;

    command_gl_ = new QGridLayout;

    perturb_pb_ = new QPushButton;
    perturb_pb_->setText(tr("Perturb !"));

    reset_pb_ = new QPushButton;
    reset_pb_->setText(tr("Reset"));

    simulate_pb_ = new QPushButton;
    simulate_pb_->setText(tr("Simulate"));

    set_input_pb_ = new QPushButton;
    set_input_pb_->setText(tr("Disturb !"));

    command_gl_->addWidget(perturb_pb_  , 0,0,1,1);
    command_gl_->addWidget(reset_pb_    , 0,1,1,1);
    command_gl_->addWidget(simulate_pb_ , 1,0,1,1);
    command_gl_->addWidget(set_input_pb_, 1,1,1,1);
    command_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Ignored,QSizePolicy::Expanding),2,2);

    command_gb_->setLayout(command_gl_);
    command_gb_->setTitle(tr("Command"));

    //-- Input
    input_gb_ = new QGroupBox;
    input_gl_ = new QGridLayout;

    voltage_in_label_ = new QLabel;
    voltage_in_label_->setText(tr("Voltage : "));
    voltage_in_dsb_ = new QDoubleSpinBox;
    voltage_in_dsb_->setSingleStep(1e-3);
    voltage_in_dsb_->setDecimals(3);
    voltage_in_dsb_->setMinimum(-24.0);
    voltage_in_dsb_->setMaximum(24.0); // 24 volt
    voltage_in_dsb_->setValue(.0);
    voltage_in_dsb_->adjustSize();

    torque_in_label_ = new QLabel;
    torque_in_label_->setText(tr("Torque : "));
    torque_in_dsb_ = new QDoubleSpinBox;
    torque_in_dsb_->setSingleStep(1e-3);
    torque_in_dsb_->setDecimals(3);
    torque_in_dsb_->setMinimum(.0);
    torque_in_dsb_->setMaximum(10.0);
    torque_in_dsb_->setValue(.0);
    torque_in_dsb_->adjustSize();

    input_gl_->addWidget(voltage_in_label_, 0,0,1,1);
    input_gl_->addWidget(voltage_in_dsb_,   0,1,1,1);
    input_gl_->addWidget(torque_in_label_,  1,0,1,1);
    input_gl_->addWidget(torque_in_dsb_,    1,1,1,1);
    input_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Ignored,QSizePolicy::Expanding),2,2);

    input_gb_->setLayout(input_gl_);
    input_gb_->setTitle(tr("Input Disturbance"));

    //-- Fault
    fault_gb_ = new QGroupBox;
    fault_gl_ = new QGridLayout;

    target_label_ = new QLabel;
    target_label_->setText(tr("Sensor : "));

    target_cb_ = new QComboBox;
    target_cb_->addItem(tr("Position"));
    target_cb_->addItem(tr("Velocity"));
    target_cb_->addItem(tr("Current"));

    fault_pb_ = new QPushButton;
    fault_pb_->setText(tr("Fault !"));

    bias_label_ = new QLabel;
    bias_label_->setText(tr("Bias : "));
    bias_dsb_ =  new QDoubleSpinBox;
    // Settings DSB params
    bias_dsb_->setValue(.0);
    bias_dial_ = new QDial;

    scale_label_ = new QLabel;
    scale_label_->setText(tr("Scale : "));
    scale_dsb_ =  new QDoubleSpinBox;
    // Settings DSB params    
    scale_dsb_->setMaximum(2.);
    scale_dsb_->setMinimum(.0);
    scale_dsb_->setValue(1.);
    scale_dsb_->setDecimals(3);
    scale_dsb_->setSingleStep(1e-3);
    scale_dial_ = new QDial;    
    scale_dial_->setMaximum(2000);
    scale_dial_->setMinimum(0);
    scale_dial_->setValue(1000);
    scale_dial_->setSingleStep(1);

    dead_zone_label_ = new QLabel;
    dead_zone_label_->setText(tr("Dead Zone : "));
    dead_zone_dsb_ =  new QDoubleSpinBox;
    // Settings DSB params
    dead_zone_dsb_->setValue(.0);
    dead_zone_dial_ = new QDial;

    fault_gl_->addWidget(target_label_,    0,0,1,1);
    fault_gl_->addWidget(target_cb_,       0,1,1,1);
    fault_gl_->addWidget(fault_pb_,      0,2,1,1);
    fault_gl_->addWidget(bias_label_,      1,0,1,1);
    fault_gl_->addWidget(bias_dsb_,        1,1,1,1);
    fault_gl_->addWidget(bias_dial_,       1,2,1,1);
    fault_gl_->addWidget(scale_label_,     2,0,1,1);
    fault_gl_->addWidget(scale_dsb_,       2,1,1,1);
    fault_gl_->addWidget(scale_dial_,      2,2,1,1);
    fault_gl_->addWidget(dead_zone_label_, 3,0,1,1);
    fault_gl_->addWidget(dead_zone_dsb_,   3,1,1,1);
    fault_gl_->addWidget(dead_zone_dial_,  3,2,1,1);
    fault_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),4,3);

    fault_gb_->setLayout(fault_gl_);
    fault_gb_->setTitle(tr("Fault"));

    //-- Watermark
    watermark_widget_ = new QWidget;
    watermark_gl_ = new QGridLayout;

    QImageReader image_reader("../example/gui/Logo_Universitas_Gadjah_Mada.png");
    logo_image_ = new QImage;
    *logo_image_ = image_reader.read().scaled(QSize(200,200), Qt::KeepAspectRatio);
    logo_image_label_ = new QLabel;
    logo_image_label_->setPixmap(QPixmap::fromImage(*logo_image_));

    title_label_ = new QLabel;
    title_label_->setText(tr("Control and Sensor Fault Isolation\n"
                             "   Design on DC Motor along with\n"
                             "           Remote Monitoring"));
    title_label_->setStyleSheet("font: bold italic");

    watermark_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding), 0,0);
    watermark_gl_->addWidget(logo_image_label_, 1,0,1,1);
    watermark_gl_->setAlignment(logo_image_label_, Qt::AlignCenter);
    watermark_gl_->addWidget(title_label_,      2,0,1,1);
    watermark_gl_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding), 3,1);
    watermark_widget_->setLayout(watermark_gl_);

    //-- Main
    main_widget_ = new QWidget;
    main_layout_ = new QGridLayout;
    main_layout_->addWidget(params_gb_,        0,0,2,1);
    main_layout_->addWidget(input_gb_,         0,1,1,1);
    main_layout_->addWidget(command_gb_,       1,1,1,1);
    main_layout_->addWidget(fault_gb_,         2,0,1,1);
    main_layout_->addWidget(watermark_widget_, 2,1,1,1);

//    main_layout_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),1,2);
    main_layout_->addItem(new QSpacerItem(0,0,QSizePolicy::Expanding,QSizePolicy::Expanding),3,2);
    main_widget_->setLayout(main_layout_);

    this->setCentralWidget(main_widget_);
    this->adjustSize();
    this->layout()->setSizeConstraint(QLayout::SetFixedSize);
    this->setWindowTitle(tr("Control Simulator"));

}

void MainWindow::setupControllerDialog(){
    controller_dialog_ = new jacl::ControllerDialog(this);
}

void MainWindow::setupMenus(){
    controller_act_ = new QAction(tr("Controller"), this);
    controller_act_->setStatusTip(tr("Open controller dialog"));
    connect(controller_act_, SIGNAL(triggered()), this, SLOT(openControllerDialog()));
    tools_menu_ = this->menuBar()->addMenu(tr("Tools"));
    tools_menu_->addAction(controller_act_);
}

void MainWindow::setupActions(){
    connect(perturb_pb_, SIGNAL(clicked()), this, SLOT(perturbAct()));
    connect(reset_pb_, SIGNAL(clicked()), this, SLOT(resetAct()));
    connect(simulate_pb_, SIGNAL(clicked()), this, SLOT(simulateAct()));
    connect(set_input_pb_, SIGNAL(clicked()), this, SLOT(setInputAct()));

    connect(bias_dsb_, SIGNAL(valueChanged(double)), this, SLOT(biasDialConv(double)));
    connect(scale_dsb_, SIGNAL(valueChanged(double)), this, SLOT(scaleDialConv(double)));
    connect(dead_zone_dsb_, SIGNAL(valueChanged(double)), this, SLOT(deadZoneDialConv(double)));

    connect(bias_dial_, SIGNAL(valueChanged(int)), this, SLOT(biasDSBConv(int)));
    connect(scale_dial_, SIGNAL(valueChanged(int)), this, SLOT(scaleDSBConv(int)));
    connect(dead_zone_dial_, SIGNAL(valueChanged(int)), this, SLOT(deadZoneDSBConv(int)));

    connect(controller_dialog_, SIGNAL(setRefSig()), this, SLOT(refAct()));
    connect(controller_dialog_, SIGNAL(setModeSig(int)), this, SLOT(modeAct(int)));
    qRegisterMetaType<arma::vec>("arma::vec");
    connect(this, SIGNAL(setDataMonitor(arma::vec)), controller_dialog_, SLOT(setDataMonitor(arma::vec)));
    qRegisterMetaType<arma::Col<uint8_t>>("arma::Col<uint8_t>");
    connect(this, SIGNAL(setSensorStatus(arma::Col<uint8_t>)), controller_dialog_, SLOT(setSensorStatus(arma::Col<uint8_t>)));
}

void MainWindow::perturbAct(){

//    for(int idx(iBm); idx <= iRa; idx++)
//        ss_.param(idx) = params_dsb_[idx]->value();

//    &bm, &jm, &ki, &la, &kb, &ra
    for(auto &pair:{std::pair<std::decay<jacl::PhysicalParameter>::type*, QDoubleSpinBox*>
                        {&bm,params_dsb_[iBm]},
                        {&jm,params_dsb_[iJm]},
                        {&ki,params_dsb_[iKi]},
                        {&la,params_dsb_[iLa]},
                        {&kb,params_dsb_[iKb]},
                        {&ra,params_dsb_[iRa]}})
        pair.first->perturbed = pair.second->value();

    // cobserver_plt_.updateVar();
//    sim_.setLinearStateSpace(ss_.A(), ss_.B(), ss_.C(), ss_.D());
}

void MainWindow::resetAct(){

    params_dsb_[iBm]->setValue(bm.nominal);
    params_dsb_[iJm]->setValue(jm.nominal);
    params_dsb_[iKi]->setValue(ki.nominal);
    params_dsb_[iLa]->setValue(la.nominal);
    params_dsb_[iKb]->setValue(kb.nominal);
    params_dsb_[iRa]->setValue(ra.nominal);

    perturbAct();
}

void MainWindow::simulateAct(){
    //-- Continuous
    // csys_plt_.start();
    // cobserver_plt_.start();
    // posctrl_plt_.start();
    // spdctrl_plt_.start();
    // cobserver_plt_.simulate();

    //-- Discrete
    dsys_simo_plt_.start();
    pos_dctrl_plt_.start();
    // dobserver_simo_plt_.start();
    sifd_.viewSignals();
}

void MainWindow::setInputAct(){
    arma::mat in(2, 1);
    in(0) += voltage_in_dsb_->value();
    in(1) += torque_in_dsb_->value();
    din_ += in.submat(0,0,0,0);
}

void MainWindow::biasDialConv(double _val){
    bias_dial_->setValue((int)_val);
}

void MainWindow::scaleDialConv(double _val){
    scale_dial_->setValue(_val * 1000);
}

void MainWindow::deadZoneDialConv(double _val){
    dead_zone_dial_->setValue((int)_val);
}

void MainWindow::biasDSBConv(int _val){
    bias_dsb_->setValue((double)_val);
}

void MainWindow::scaleDSBConv(int _val){
    scale_dsb_->setValue((double)_val * .001);
}

void MainWindow::deadZoneDSBConv(int _val){
    dead_zone_dsb_->setValue((double)_val);
}

void MainWindow::openControllerDialog(){
    controller_dialog_->show();
}

void MainWindow::refAct(){
    ref_(0) = controller_dialog_->getPosition();
    ref_(1) = controller_dialog_->getVelocity();
    ref_(2) = controller_dialog_->getCurrent();
}

void MainWindow::modeAct(int _val){
    control_mode_ = _val;
}

void MainWindow::makeFault(arma::vec* _out){
    int target_idx = target_cb_->currentIndex();
    (*_out)(target_idx) += bias_dsb_->value();
    if(std::fabs((*_out)(target_idx)) < dead_zone_dsb_->value())
        (*_out)(target_idx) = 0;
    (*_out)(target_idx) *= scale_dsb_->value();
}

void MainWindow::setupSIMODCMotor(){
    //-- SIMO DC motor open-loop
    {
        SIMO::Formulas fA{
            JC(simo_,.0),JC(simo_,1.0),JC(simo_,.0),
            JC(simo_,.0),JE(simo_,-_p(iBm)/_p(iJm)),JE(simo_,_p(iKi)/_p(iJm)),
            JC(simo_,.0),JE(simo_,-_p(iKb)/_p(iLa)),JE(simo_,-_p(iRa)/_p(iLa))
        };

        SIMO::Formulas fB{
            JC(simo_, .0),
            JC(simo_, .0),
            JE(simo_, 1.0/_p(iLa))
        };

        SIMO::Formulas fC{
            JC(simo_, 1.0), JC(simo_, .0), JC(simo_, .0),
            JC(simo_, .0), JC(simo_, 1.0), JC(simo_, .0),
            JC(simo_, .0), JC(simo_, .0), JC(simo_, 1.0)
        };

        SIMO::Formulas fD{
            JC(simo_, .0),
            JC(simo_, .0),
            JC(simo_, .0)
        };
        simo_.setA(fA);
        simo_.setB(fB);
        simo_.setC(fC);
        simo_.setD(fD);
    }
    jacl::parser::saveStateSpace(simo_, "motor_dc_simo.jacl");
    jacl::common::StateSpacePack dsimo = jacl::common::discretize(simo_, SAMPLING_PERIOD);
    dsimo_.setA(std::get<0>(dsimo)); dsimo_.setB(std::get<1>(dsimo));
    dsimo_.setC(std::get<2>(dsimo)); dsimo_.setD(std::get<3>(dsimo));
    dsimo_.A().print("Ad : "); dsimo_.B().print("Bd : ");
    dsimo_.C().print("Cd : "); dsimo_.D().print("Dd : ");

    arma::mat dobsv_gain;
    arma::mat dpoles{-0.8,-0.85,0.83};
    jacl::pole_placement::KautskyNichols(&dsimo_, dpoles, &dobsv_gain, jacl::pole_placement::PolePlacementType::Observer);

    dobserver_simo_.setGain(dobsv_gain.t());
    dobsv_gain.print("\nDiscrete Observer Gain : ");
    jacl::parser::saveGain(dobsv_gain, "discrete_observer_gain.jacl");
//    jacl::parser::readGain(&obsv_gain, "discrete_observer_gain.jacl");

    dsys_simo_plt_.init();
    dsys_simo_plt_.setTitle("SIMO DC motor");
    dsys_simo_plt_.setPlotName({"Position", "Velocity", "Current", "Voltage"});

    dobserver_simo_plt_.init();
    dobserver_simo_plt_.setTitle("Observer SIMO DC motor");
    dobserver_simo_plt_.setPlotName({"Est. Position", "Est. Velocity", "Est. Current"});

    {
        MReal::Formulas fA{
            //-- Row-1
            JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-2
            JC(m_real_, .0), JE(m_real_, -_p(iBm)/_p(iJm)), JE(m_real_, _p(iKi)/_p(iJm)),
            JE(m_real_, weight_(iBm)*-1./_p(iJm)), JE(m_real_, weight_(iJm)*-1./_p(iJm)), JE(m_real_, weight_(iKi)*-1./_p(iJm)), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            // JE(m_real_, -1./_p(iJm)), JE(m_real_, -1./_p(iJm)), JE(m_real_, -1./_p(iJm)), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-3
            JC(m_real_, .0), JE(m_real_, -_p(iKb)/_p(iLa)), JE(m_real_, -_p(iRa)/_p(iLa)),
            JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JE(m_real_, weight_(iLa)*-1./_p(iLa)), JE(m_real_, weight_(iKb)*-1./_p(iLa)), JE(m_real_, weight_(iRa)*-1./_p(iLa)),
            // JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JE(m_real_, -1./_p(iLa)), JE(m_real_, -1./_p(iLa)), JE(m_real_, -1./_p(iLa)),
            //-- Row-4
            JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-5
            JC(m_real_, .0), JE(m_real_, -_p(iBm)/_p(iJm)), JE(m_real_, _p(iKi)/_p(iJm)),
            JE(m_real_, weight_(iBm)*-1./_p(iJm)), JE(m_real_, weight_(iJm)*-1./_p(iJm)), JE(m_real_, weight_(iKi)*-1./_p(iJm)), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            // JE(m_real_, -1./_p(iJm)), JE(m_real_, -1./_p(iJm)), JE(m_real_, -1./_p(iJm)), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-6
            JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-7
            JC(m_real_, .0), JE(m_real_, -_p(iKb)/_p(iLa)), JE(m_real_, -_p(iRa)/_p(iLa)),
            JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JE(m_real_, weight_(iLa)*-1./_p(iLa)), JE(m_real_, weight_(iKb)*-1./_p(iLa)), JE(m_real_, weight_(iRa)*-1./_p(iLa)),
            // JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JE(m_real_, -1./_p(iLa)), JE(m_real_, -1./_p(iLa)), JE(m_real_, -1./_p(iLa)),
            //-- Row-8
            JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0),
            //-- Row-9
            JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0)
        };

        MReal::Formulas fB{
            JC(m_real_, .0),        
            JC(m_real_, .0),
            JE(m_real_, -1./_p(iLa)),
            JC(m_real_, .0),
            JC(m_real_, .0),
            JC(m_real_, .0),
            JE(m_real_, -1./_p(iLa)),
            JC(m_real_, .0),
            JC(m_real_, .0)
        };

        MReal::Formulas fC{
            JC(m_real_, 1.), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, .0)
            // JC(m_real_, .0), JC(m_real_, 1.), JC(m_real_, .0),
            // JC(m_real_, .0), JC(m_real_, .0), JC(m_real_, 1.)
        };

        MReal::Formulas fD{
            JC(m_real_, .0)
            // JC(m_real_, .0),
            // JC(m_real_, .0)
        };
        m_real_.setA(fA);
        m_real_.setB(fB);
        m_real_.setC(fC);
        m_real_.setD(fD);
    }
    jacl::common::StateSpacePack dm_real_ = jacl::common::discretize(m_real_, SAMPLING_PERIOD);
    {
        arma::mat temp1, temp2, temp3;
        arma::mat zeros9x1(9,1,arma::fill::zeros);
        arma::mat zeros1x9(1,9,arma::fill::zeros);
        arma::mat zeros1x1(1,1,arma::fill::zeros);
        arma::mat eye1x1(1,1,arma::fill::eye);

        temp1 = arma::join_horiz(zeros9x1, std::get<1>(dm_real_));
        temp2 = arma::join_horiz(zeros9x1, temp1);
        arma::mat B = arma::join_horiz(std::get<1>(dm_real_), temp2);

        temp1 = arma::join_vert(zeros1x9, -std::get<2>(dm_real_));
        arma::mat C = arma::join_vert(-std::get<2>(dm_real_), temp1);

        temp1 = arma::join_horiz(zeros1x1, eye1x1);
        temp2 = arma::join_horiz(zeros1x1, temp1);
        temp3 = arma::join_horiz(zeros1x1, arma::join_horiz(zeros1x1, zeros1x1));
        arma::mat D11 = arma::join_vert(temp2, temp3);

        arma::mat D12 = arma::join_vert(zeros1x1, eye1x1);

        temp1 = arma::join_horiz(-eye1x1,eye1x1);
        arma::mat D21 = arma::join_horiz(zeros1x1,temp1);
        
        arma::mat D22 = zeros1x1;

        temp1 = arma::join_horiz(D11, D12);
        temp2 = arma::join_horiz(D21, D22);
        arma::mat D = arma::join_vert(temp1, temp2);

        m_icm_.setA(std::get<0>(dm_real_));
        m_icm_.setB(B);
        m_icm_.setC(C);
        m_icm_.setD(D);
    }

    {
        pos_real_.setA(dsimo_.A());
        pos_real_.setB(dsimo_.B());
        pos_real_.setC({1.,0.,0.});
        pos_real_.setD(arma::zeros(1,1));

        arma::mat temp1, temp2, temp3;
        arma::mat zeros3x1(3,1,arma::fill::zeros);
        arma::mat zeros1x3(1,3,arma::fill::zeros);
        arma::mat zeros1x1(1,1,arma::fill::zeros);
        arma::mat eye1x1(1,1,arma::fill::eye);

        temp1 = arma::join_horiz(zeros3x1, pos_real_.B());
        temp2 = arma::join_horiz(zeros3x1, temp1);
        arma::mat B = arma::join_horiz(pos_real_.B(), temp2);

        temp1 = arma::join_vert(zeros1x3, -pos_real_.C());
        arma::mat C = arma::join_vert(-pos_real_.C(), temp1);

        temp1 = arma::join_horiz(zeros1x1, eye1x1);
        temp2 = arma::join_horiz(zeros1x1, temp1);
        temp3 = arma::join_horiz(zeros1x1, arma::join_horiz(zeros1x1, zeros1x1));
        arma::mat D11 = arma::join_vert(temp2, temp3);

        arma::mat D12 = arma::join_vert(zeros1x1, eye1x1);

        temp1 = arma::join_horiz(-eye1x1,eye1x1);
        arma::mat D21 = arma::join_horiz(zeros1x1,temp1);
        
        arma::mat D22 = zeros1x1;

        temp1 = arma::join_horiz(D11, D12);
        temp2 = arma::join_horiz(D21, D22);
        arma::mat D = arma::join_vert(temp1, temp2);

        pos_icm_.setA(pos_real_.A());
        pos_icm_.setB(B);
        pos_icm_.setC(C);
        pos_icm_.setD(D);
    }
    pos_icm_.A().print("Aicm : ");
    pos_icm_.B().print("Bicm : ");
    pos_icm_.C().print("Cicm : ");
    pos_icm_.D().print("Dicm : ");
    pos_dhinf_ = new PosDHinf(&pos_sys_, 1.7);
    auto K( pos_dhinf_->solve() );
    pos_dctrl_.setA(std::get<0>(K)); pos_dctrl_.setB(std::get<1>(K));
    pos_dctrl_.setC(std::get<2>(K)); pos_dctrl_.setD(std::get<3>(K));

    //-- Closed-loop
    jacl::state_space::Linear<double, 6,1,1> cl_ss;
    jacl::system::Discrete<jacl::state_space::Linear<double, 6,1,1>> cl_sys(&cl_ss);
    {
        arma::mat temp1, temp2;

        temp1 = pos_real_.B()*pos_dctrl_.D()*pos_real_.C();
        arma::mat A = arma::join_cols(
            arma::join_rows(pos_real_.A() - temp1, pos_real_.B()*pos_dctrl_.C()),
            arma::join_rows(-pos_dctrl_.B()*pos_real_.C(), pos_dctrl_.A())
        );
        
        arma::mat B = arma::join_cols(
            pos_real_.B()*pos_dctrl_.D(),
            pos_dctrl_.B()
        );

        cl_ss.setA(A);
        cl_ss.setB(B);
        cl_ss.setC(arma::join_rows(pos_real_.C(), arma::zeros(1,3)));
        cl_ss.setD(pos_real_.D());
    }
    cl_ss.A().print("Acl : ");
    cl_ss.B().print("Bcl : ");
    cl_ss.C().print("Ccl : ");
    cl_ss.D().print("Dcl : ");
    jacl::analysis::transient_data_t transient_data = jacl::analysis::transient(cl_sys);
    std::cout << "Rise time : " << jacl::analysis::getRiseTime(transient_data) << std::endl;
    std::cout << "Peak time : " << jacl::analysis::getPeakTime(transient_data) << std::endl;
    std::cout << "Overshoot : " << jacl::analysis::getOvershoot(transient_data) << std::endl;
    std::cout << "Settling time : " << jacl::analysis::getSettlingTime(transient_data) << std::endl;

    pos_dctrl_plt_.init();
    pos_dctrl_plt_.setTitle("Discrete Position Controller");
    pos_dctrl_plt_.setPlotName({"Error","Input"});
}

void MainWindow::setupPositionController(){

    siso_pos_.setA({{0.,1.,0.},
                    {0.,-.6667,722.2},
                    {0,-83.87,-3339}});
    siso_pos_.setB(arma::colvec({0.,
                    0.,
                    1613.}));
    siso_pos_.setC({1.,0.,0.});
    siso_pos_.setD(arma::colvec({0.}));

    InterConnMatPos::Formulas fA = {
        JC(icm_pos_, -1.444e+04), JC(icm_pos_, -7780), JC(icm_pos_, -2912), JC(icm_pos_, -841.9), JC(icm_pos_,-269), JC(icm_pos_, 0),
        JC(icm_pos_, 8192), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 4096), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 1024), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 256), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 1), JC(icm_pos_, 0)
    };


    InterConnMatPos::Formulas fB = {
        JC(icm_pos_, 64), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 64),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0)
    };

    InterConnMatPos::Formulas fC = {
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 2.276e-06), JC(icm_pos_, 0.04105), JC(icm_pos_, 0.5249), JC(icm_pos_, 77.95),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 2.276e-06), JC(icm_pos_, 0.04105), JC(icm_pos_, 0.5249), JC(icm_pos_, 77.95)
    };

    InterConnMatPos::Formulas fD = {
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 1), JC(icm_pos_, 0),
        JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 0), JC(icm_pos_, 1),
        JC(icm_pos_, 0), JC(icm_pos_, -1), JC(icm_pos_, 1), JC(icm_pos_, 0)
    };    

    icm_pos_.setA(fA);
    icm_pos_.setB(fB);
    icm_pos_.setC(fC);
    icm_pos_.setD(fD);
    jacl::system::Continuous<InterConnMatPos> icm_pos_sys(&icm_pos_);
    hinf_pc_ = new HInfPC(&icm_pos_sys, 2.7882);
    jacl::common::StateSpacePack K( hinf_pc_->solve() );
    k_pos_.setA(std::get<0>(K));
    k_pos_.setB(std::get<1>(K));
    k_pos_.setC(std::get<2>(K));
    k_pos_.setD(std::get<3>(K));    

    k_pos_.A().print("A : ");
    k_pos_.B().print("B : ");
    k_pos_.C().print("C : ");
    k_pos_.D().print("D : ");

    jacl::parser::saveStateSpace(k_pos_, "position_controller.jacl");

    //-- Example read
    /*PosCtrl controller;
    jacl::parser::readStateSpace(&controller, "position_controller.jacl");
    controller.A().print("Controller A : ");
    controller.B().print("Controller B : ");
    controller.C().print("Controller C : ");
    controller.D().print("Controller D : ");*/

    posctrl_plt_.init();
    posctrl_plt_.setTitle("Position Controller");
    posctrl_plt_.setDelay() = SAMPLING_PERIOD;
    posctrl_plt_.setPlotName({"Position Error", "Voltage"});
}

void MainWindow::setupSpeedController(){
    InterConnMatSpd::Formulas fA = {
        JC(icm_spd_, -1.444e+04), JC(icm_spd_, -7780), JC(icm_spd_, -2912), JC(icm_spd_, -841.9), JC(icm_spd_, -269),
        JC(icm_spd_, 8192), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 4096), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 1024), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 256), JC(icm_spd_, 0)
    };

    InterConnMatSpd::Formulas fB = {
        JC(icm_spd_, 128), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 128),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0)
    };

    InterConnMatSpd::Formulas fC = {
        JC(icm_spd_, 0), JC(icm_spd_, 0.004661), JC(icm_spd_, 24.99), JC(icm_spd_, 80.17), JC(icm_spd_, 48.49),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0.004661), JC(icm_spd_, 24.99), JC(icm_spd_, 80.17), JC(icm_spd_, 48.49)
    };

    InterConnMatSpd::Formulas fD = {
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 1), JC(icm_spd_, 0),
        JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 0), JC(icm_spd_, 1),
        JC(icm_spd_, 0), JC(icm_spd_, -1), JC(icm_spd_, 1), JC(icm_spd_, 0)
    };

    icm_spd_.setA(fA);
    icm_spd_.setB(fB);
    icm_spd_.setC(fC);
    icm_spd_.setD(fD);
    jacl::system::Continuous<InterConnMatSpd> icm_spd_sys(&icm_spd_);
    hinf_sc_ = new HInfSC(&icm_spd_sys, 3.1);
    jacl::common::StateSpacePack K( hinf_sc_->solve() );
    k_spd_.setA(std::get<0>(K));
    k_spd_.setB(std::get<1>(K));
    k_spd_.setC(std::get<2>(K));
    k_spd_.setD(std::get<3>(K));

    jacl::parser::saveStateSpace(k_spd_, "speed_controller.jacl");

    k_spd_.A().print("Kspd_A : ");
    k_spd_.B().print("Kspd_B : ");
    k_spd_.C().print("Kspd_C : ");
    k_spd_.D().print("Kspd_D : ");

    spdctrl_plt_.init();
    spdctrl_plt_.setTitle("Speed Controller");
    spdctrl_plt_.setDelay() = SAMPLING_PERIOD;
    spdctrl_plt_.setPlotName({"Velocity Error", "Voltage"});
}

void MainWindow::setupNLP(){
    {        
        NLPendulum::Formulas state = {
            JE(nlp_, _p(ix2)),
            JE(nlp_, -1*(_p(ipg)/_p(ipl))*sin(_p(ix1)) - (_p(ipk)/_p(ipm))*_p(ix1)), 
        };
        NLPendulum::Formulas output = {
            JE(nlp_, _p(ix1))
        };
    }
}

double MainWindow::angularSpeed2Voltage(double _speed, double _torque){
    return _torque*ra.perturbed/ki.perturbed + ki.perturbed*_speed;
}

void MainWindow::closedLoopProcess(){
    //-- Continuous
    arma::vec out(3,1,arma::fill::zeros),in(1,1,arma::fill::zeros);
    arma::vec err(in),est(out);
    din_ = in;
    control_mode_ = jacl::traits::toUType(jacl::ControllerDialog::ControlMode::Position);
//    arma::mat last_err(err);
//    arma::mat diff(err);
//    auto Kp(10.), Kd(1.);
    // ifd_.init({{-.76,-.65}, {-.63,-.51}, {-.86,-.72}});
    sifd_.init({{-.11,-.15}}, "SIFD", {"Est. Curr.","Est. Spd.","Est. Pos."});
    arma::cx_vec p = jacl::common::poles(dsimo_);    
    p.print("Discrete DC motor poles : ");
    p = jacl::common::poles(simo_);
    p.print("Continuous DC motor poles : ");        
    while(cl_status_){
        //-- PD Control
       /*err(1) = ref_(1) - out(1);
       diff(1) = err(1) - last_err(1);
       last_err(1) = err(1);
       err.print("Error : ");
       in(0) = angularSpeed2Voltage(Kp*err(1) + Kd*diff(1),.0);
       if(in(0) > 12.0){
           in(0) = 12.0;
       }
       if(in(0) < -12.0){
           in(0) = -12.0;
       }*/

        //-- Error between reference and output
        err(control_mode_) = ref_(control_mode_) - out(control_mode_);

        //-- H-infinity Controller
        // if(control_mode_ == jacl::traits::toUType(jacl::ControllerDialog::ControlMode::Position)){
        //     err(1) = 0;
        //     in.submat(0,0,0,0) = posctrl_sys_.convolve(err.submat(0,0,0,0));
        // }else if(control_mode_ == jacl::traits::toUType(jacl::ControllerDialog::ControlMode::Velocity)){
        //     err(0) = 0;
        //     in.submat(0,0,0,0) = spdctrl_sys_.convolve(err.submat(1,0,1,0));
        // }else{
        //     //-- do nothing
        // }

        //-- Continuous
        // out = csys_.convolve(in);
        // est = cobserver_.convolve(in, out);

        //-- Discrete
        din_ = pos_dctrl_sys_.convolve(err);

        //-- saturation
        if(din_(0) > 24.)
            din_(0) = 24.;
        else if(din_(0) < -24.)
            din_(0) = -24.;
        
        out = dsys_simo_.convolve(din_);
        makeFault(&out);
        SIFD::diag_pack_t diag_pack = sifd_.detect(din_, out);        
        Q_EMIT setDataMonitor(
            arma::join_cols(out,
                arma::join_cols(
                    std::get<0>(diag_pack[0]),
                    arma::join_cols(std::get<0>(diag_pack[1]),
                                    std::get<0>(diag_pack[2])
                    )
                )
            )
        );
        arma::Col<uint8_t> status{std::get<1>(diag_pack[0]), std::get<1>(diag_pack[1]), std::get<1>(diag_pack[2])};
        Q_EMIT setSensorStatus(status);
        boost::this_thread::sleep_for(boost::chrono::milliseconds((int)(SAMPLING_PERIOD*1000.)));
    }
}
