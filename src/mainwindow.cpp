/**
*   @author : koseng (Lintang)
*   @brief : GUI Source
*/

#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
    , bm(.000048)
    , jm(7.2e-6)
    , ki(.052)
    , la(.00062)
    , kb(.052)
    , ra(2.07)
    , ss_(&bm, &jm, &ki, &la, &kb, &ra)
    , G_(&bm, &jm, &ki, &la, &kb, &ra)
//    , sim_(3,2,2)
    , system_sim_(&ss_)
    , obs_gain_({{ .0, .0, .0},
                 { .0, .0, .0},
                 { .0, .0, .0}})
    , observer_sim_(&ss_, obs_gain_)
    , controller_sim_(&K_){

    ui->setupUi(this);

    this->removeToolBar(ui->mainToolBar);

    QFile style_file("../gui/dark_style_lintang.qss");
    style_file.open(QFile::ReadOnly);
    QString style_sheet(style_file.readAll());
    this->setStyleSheet(style_sheet);

    //-- DC Motor Open-Loop
    std::cout << "Preparing system ..." << std::endl;
    {
        StateSpace::Formula A11 = JC(ss_,.0); StateSpace::Formula A12 = JC(ss_,1.0);
        StateSpace::Formula A13 = JC(ss_,.0); StateSpace::Formula A21 = JC(ss_,.0);
        StateSpace::Formula A22 = JE(ss_,-ss_(iBm)/ss_(iJm)); StateSpace::Formula A23 = JE(ss_,ss_(iKi)/ss_(iJm));
        StateSpace::Formula A31 = JC(ss_,.0); StateSpace::Formula A32 = JE(ss_,-ss_(iKb)/ss_(iLa));
        StateSpace::Formula A33 = JE(ss_,-ss_(iRa)/ss_(iLa));

        StateSpace::Formulas fA{
            A11,A12,A13,
            A21,A22,A23,
            A31,A32,A33
        };

        StateSpace::Formula B11 = JC(ss_, .0); StateSpace::Formula B12 = JC(ss_, .0);
        StateSpace::Formula B21 = JC(ss_, .0); StateSpace::Formula B22 = JE(ss_, -1.0/ss_(iJm));
        StateSpace::Formula B31 = JE(ss_, 1.0/ss_(iLa)); StateSpace::Formula B32 = JC(ss_, .0);

        StateSpace::Formulas fB{
            B11, B12,
            B21, B22,
            B31, B32
        };

        StateSpace::Formula C11 = JC(ss_, 1.0); StateSpace::Formula C12 = JC(ss_, .0);
        StateSpace::Formula C13 = JC(ss_, .0); StateSpace::Formula C21 = JC(ss_, .0);
        StateSpace::Formula C22 = JC(ss_, 1.0); StateSpace::Formula C23 = JC(ss_, .0);
        StateSpace::Formula C31 = JC(ss_, .0); StateSpace::Formula C32 = JC(ss_, .0);
        StateSpace::Formula C33 = JC(ss_, 1.0);

        StateSpace::Formulas fC{
            C11, C12, C13,
            C21, C22, C23,
            C31, C32, C33
        };

        StateSpace::Formula D11 = JC(ss_, .0); StateSpace::Formula D12 = JC(ss_, .0);
        StateSpace::Formula D21 = JC(ss_, .0); StateSpace::Formula D22 = JC(ss_, .0);
        StateSpace::Formula D31 = JC(ss_, .0); StateSpace::Formula D32 = JC(ss_, .0);

        StateSpace::Formulas fD{
            D11, D12,
            D21, D22,
            D31, D32,
        };

        ss_.setA(fA);
        ss_.setB(fB);
        ss_.setC(fC);
        ss_.setD(fD);
        ss_.formulaToMat();
    }

    ss_.A().print("A : ");
    ss_.B().print("B : ");
    ss_.C().print("C : ");
    ss_.D().print("D : ");

    std::cout << "Ctrb : " << jacl::common::controllable(ss_.A(), ss_.B()) << std::endl;
    std::cout << "Obsv : " << jacl::common::observable(ss_.A(), ss_.C()) << std::endl;


    //-- G Realization
    {
        GRealization::Formula A11, A12, A13;
        A11 = JC(G_, .0); A12 = JC(G_, 1.0); A13 = JC(G_, .0);
        GRealization::Formula A21, A22, A23;
        A21 = JC(G_, .0); A22 = JE(G_, -G_(iBm)/G_(iJm)); A23 = JE(G_, G_(iKi)/G_(iJm));
        GRealization::Formula A31, A32, A33;
        A31 = JC(G_, .0); A32 = JE(G_, -G_(iKb)/G_(iLa)); A33 = JE(G_, -G_(iRa)/G_(iLa));

        GRealization::Formula B11, B12, B13, B14, B15, B16, B17, B18;
        B11 = B12 = B13 = B14 = B15 = B16 = B17 = B18 = JC(G_, .0);
        GRealization::Formula B21, B22, B23, B24, B25, B26, B27, B28;
        B21 = B22 = B23 = JE(G_, -1.0/G_(iJm)); B24 = B25 = B26 = B27 = JC(G_, .0); B28 = JE(G_, -1.0/G_(iJm));
        GRealization::Formula B31, B32, B33, B34, B35, B36, B37, B38;
        B31 = B32 = B33 = JC(G_, .0); B34 = B35 = B36 = B37 = JE(G_, -1.0/G_(iLa)); B38 = JC(G_, .0);

        GRealization::Formula C11, C12, C13;
        C11 = JC(G_, .0); C12 = JC(G_, 1.0); C13 = JC(G_, .0);
        GRealization::Formula C21, C22, C23;
        C21 = JC(G_, .0); C22 = JE(G_, -G_(iBm)/G_(iJm)); C23 = JE(G_, G_(iKi)/G_(iJm));
        GRealization::Formula C31, C32, C33;
        C31 = JC(G_, .0); C32 = JC(G_, .0); C33 = JC(G_, 1.0);
        GRealization::Formula C41, C42, C43;
        C41 = JC(G_, .0); C42 = JE(G_, -G_(iKb)/G_(iLa)); C43 = JE(G_, -G_(iRa)/G_(iLa));
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
        D21 = D22 = D23 = JE(G_, -1/G_(iJm)); D24 = D25 = D26 = JC(G_, .0); D27 = JE(G_, -1/G_(iJm)); D28 = JC(G_, .0);
        GRealization::Formula D31, D32, D33, D34, D35, D36, D37, D38;
        D31 = D32 = D33 = D34 = D35 = D36 = D37 = D38 = JC(G_, .0);
        GRealization::Formula D41, D42, D43, D44, D45, D46, D47, D48;
        D41 = D42 = D43 = JC(G_, .0); D44 = D45 = D46 = JE(G_, -1/G_(iLa)); D47 = JC(G_, .0); D48 = JE(G_, -1/G_(iLa));
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
        G_.formulaToMat();
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
        P_.formulaToMat();
    }

    //-- Realization of Inter-connection Matrix    
    {
        arma::mat temp1, temp2, temp3;
        arma::mat zeros9x3(9,3,arma::fill::zeros);
        arma::mat zeros9x2(9,2,arma::fill::zeros);
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

    h_inf_ = new HInf(&ICM_);
    std::tuple<arma::mat, arma::mat, arma::mat, arma::mat> K( h_inf_->solve() );
    K_.setA(std::get<0>(K));
    K_.setB(std::get<1>(K));
    K_.setC(std::get<2>(K));
    K_.setD(std::get<3>(K));

    K_.A().print("K_A : ");
    K_.B().print("K_B : ");
    K_.C().print("K_C : ");
    K_.D().print("K_D : ");

    controller_sim_.init();
    controller_sim_.setTitle("Controller");
    controller_sim_.setDelay() = TIME_STEP;
    controller_sim_.setPlotName({"1","2","3","4","5"});
    controller_sim_.updateVariables();

    //-- Test KautskyNichols
    arma::mat observer_K;
    arma::mat poles{-50,-51,-52};
    jacl::pole_placement::KautskyNichols(&ss_, poles, &observer_K, jacl::pole_placement::Observer);

    observer_K.print("\nObserver Gain : ");

    observer_sim_.setGain(observer_K.t());

    system_sim_.init();
    system_sim_.setTitle("DC Motor");
    system_sim_.setDelay() = TIME_STEP;
    system_sim_.setPlotName({"Angular Position", "Angular Velocity", "Current",
                       "Voltage In", "Torque In",
                       "Angular Position", "Angular Velocity", "Current"});
    system_sim_.updateVariables();

    observer_sim_.init();
    observer_sim_.setTitle("Full-order Luenberger Observer of DC Motor");
    observer_sim_.setDelay() = TIME_STEP;
    observer_sim_.setPlotName({"Est. Position", "Est. Velocity", "Est. Current"
                              ,"Est. Out Position", "Est. Out Velocity", "Est. Out Current"});
    observer_sim_.updateVariables();

//    sim_.init();
//    sim_.setTitle("DC Motor Simulation");
//    sim_.setDelay() = .0;
//    sim_.setPlotName({"Angular Position", "Angular Velocity", "Current",
//                      "Torque In", "Voltage In",
//                      "Angular Position", "Angular Velocity"});
//    sim_.setStateSpace(ss_.A(), ss_.B(), ss_.C(), ss_.D());    

    setupWidgets();
    setupControllerDialog();
    setupMenus();
    setupActions();

    cl_status_ = true;
    cl_thread_ = boost::thread{boost::bind(&MainWindow::closedLoopProcess, this)};

    jacl::StateSpace<12,4,4> another_ss;

}

MainWindow::~MainWindow(){
    boost::mutex::scoped_lock lk(cl_mtx_);
    cl_status_ = true;
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
    set_input_pb_->setText(tr("Set Input"));

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
    input_gb_->setTitle(tr("Input"));

    //-- Fault

    fault_gb_ = new QGroupBox;

    fault_gl_ = new QGridLayout;

    target_label_ = new QLabel;
    target_label_->setText(tr("Sensor : "));

    target_cb_ = new QComboBox;
    target_cb_->addItem(tr("Position"));
    target_cb_->addItem(tr("Velocity"));
    target_cb_->addItem(tr("Current"));

    details_pb_ = new QPushButton;
    details_pb_->setText(tr("Details"));

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
    scale_dsb_->setValue(.0);
    scale_dial_ = new QDial;

    dead_zone_label_ = new QLabel;
    dead_zone_label_->setText(tr("Dead Zone : "));
    dead_zone_dsb_ =  new QDoubleSpinBox;
    // Settings DSB params
    dead_zone_dsb_->setValue(.0);
    dead_zone_dial_ = new QDial;

    fault_gl_->addWidget(target_label_,    0,0,1,1);
    fault_gl_->addWidget(target_cb_,       0,1,1,1);
    fault_gl_->addWidget(details_pb_,      0,2,1,1);
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

    QImageReader image_reader("../gui/Logo_Universitas_Gadjah_Mada.png");
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

    connect(controller_dialog_, SIGNAL(setRefSig()), this, SLOT())
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


    ss_.formulaToMat();
    system_sim_.updateVariables();
    observer_sim_.updateVariables();
//    sim_.setStateSpace(ss_.A(), ss_.B(), ss_.C(), ss_.D());
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

//    sim_.simulate();
    system_sim_.simulate();
    controller_sim_.simulate();
    observer_sim_.simulate();
}

void MainWindow::setInputAct(){
    arma::mat in(2, 1);
    in(0) = voltage_in_dsb_->value();
    in(1) = torque_in_dsb_->value();
//    sim_.setInput(in);
    system_sim_.setInput(in);
    observer_sim_.setInput(in);
}

void MainWindow::biasDialConv(double _val){
    bias_dial_->setValue((int)_val);
}

void MainWindow::scaleDialConv(double _val){
    scale_dial_->setValue((int)_val);
}

void MainWindow::deadZoneDialConv(double _val){
    dead_zone_dial_->setValue((int)_val);
}

void MainWindow::biasDSBConv(int _val){
    bias_dsb_->setValue((double)_val);
}

void MainWindow::scaleDSBConv(int _val){
    scale_dsb_->setValue((double)_val);
}

void MainWindow::deadZoneDSBConv(int _val){
    dead_zone_dsb_->setValue((double)_val);
}

void MainWindow::openControllerDialog(){
    controller_dialog_->show();
}

void MainWindow::closedLoopProcess(){
    arma::mat a(3,1,arma::fill::zeros),b(2,1,arma::fill::zeros);
    while(cl_status_){
        a.print("a : ");
        controller_sim_.setInput(a);
        b = controller_sim_.getOutputSig();
        b.print("b : ");
        system_sim_.setInput(b);
        a = system_sim_.getOutputSig();
        a.print("a' : ");
        boost::this_thread::sleep_for(boost::chrono::milliseconds((int)TIME_STEP*1000));
    }
}
