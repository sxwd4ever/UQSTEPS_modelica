within Steps.Components;

model MockPCHeatExchangerV2
  "Mock Printed Circuit based Heat Exchanger - version 2"
  extends BaseExchanger;
  
  //extends ExchangerContainer;  

  import CP = Steps.Utilities.CoolProp; 
  import TB = Modelica.Blocks.Tables;  
  import UTIL = Modelica.Utilities;
  import MyUtil = Steps.Utilities.Util;
  
  replaceable package PBMedia = Steps.Media.SCO2; 
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_4a_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4a_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4a - column a in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  Modelica.Blocks.Types.ExternalCombiTable1D table_4a_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4a_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4a - column b in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_4b_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4b_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4b - column a in Kim[2012] for pitch=24.6, dh=1.222 (dc=1.3 mm))";  
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_4b_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4b_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4b - column b default table";  
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_4c_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4c_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4c - column a in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  Modelica.Blocks.Types.ExternalCombiTable1D table_4c_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4c_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 4c - column b in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";  

  Modelica.Blocks.Types.ExternalCombiTable1D table_5a_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5a_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5a - column c in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  Modelica.Blocks.Types.ExternalCombiTable1D table_5a_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5a_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5a - column d in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_5b_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5b_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5b - column c in Kim[2012] for pitch=24.6, dh=1.222 (dc=1.3 mm))";
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_5b_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5b_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5b - column d default table";
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_5c_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5c_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5c - column c in Kim[2012] for pitch=24.6, dh=0.922 (dc=1.3 mm))";

  Modelica.Blocks.Types.ExternalCombiTable1D table_5c_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5c_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "Table 5c - column d in Kim[2012] for pitch=12.3, dh=0.922 (dc=1.3 mm))";
    
  Modelica.Blocks.Types.ExternalCombiTable1D table_a = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4d_a", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_b = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "4d_b", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_c = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5d_c", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  
  Modelica.Blocks.Types.ExternalCombiTable1D table_d = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "5d_d", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/kim_2012.txt"), table = fill(0.0, 9, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2);
  
  inner Modelica.Blocks.Types.ExternalCombiTable1D table_th_inconel_750 = Modelica.Blocks.Types.ExternalCombiTable1D(tableName = "inconel_750", fileName = Modelica.Utilities.Files.loadResource("modelica://Steps/Resources/Data/th_conductivity.txt"), table = fill(0.0, 6, 2), smoothness = Modelica.Blocks.Types.Smoothness.LinearSegments, columns = 2:2) "thermal conductivity for inconel_750";	
  
  parameter Integer N_seg = 1 "Number of segments in a tube";
  
  parameter Modelica.SIunits.Angle phi = 0.0 "unit rad";
  
  parameter Modelica.SIunits.ReynoldsNumber Re_design = 0.0 "On-design ReynoldsNumber";
  
  parameter Modelica.SIunits.Diameter d_c = 0.0 "Diameter of semi-circular channel";
  // inlet/outlet temperature of hot/cool channel
  parameter Modelica.SIunits.Temp_C T_hot_in = 500;
  //parameter Modelica.SIunits.Temp_C T_hot_out;
  parameter Modelica.SIunits.Temp_C T_cool_in = 300;
  //parameter Modelica.SIunits.Temp_C T_cool_out;
  
  //pressure of hot/cool channel
  parameter Modelica.SIunits.Pressure p_hot = 20 * 1e6;
  parameter Modelica.SIunits.Pressure p_cool = 9 * 1e6;
  
  // mass flow of hot/cool channel
  parameter Modelica.SIunits.MassFlowRate m_dot_hot = 5;
  parameter Modelica.SIunits.MassFlowRate m_dot_cool = 5;
  
  parameter Modelica.SIunits.Length pitch = 10 "pitch length of channel";
  
  inner parameter String name_material = "inconel 750";
  
  parameter Boolean debug_mode = false;  
  
  inner Modelica.SIunits.Diameter d_h = 4 * A_c / peri_c "Hydraulic Diameter";
  
  inner Modelica.SIunits.Length peri_c = d_c * Modelica.Constants.pi /2 + d_c "perimeter of semi-circular";
  
  inner Modelica.SIunits.Length t_wall "thickness of wall between two neighboring hot and cold";
  
  inner Integer N_channel "number of channels";
  
  inner Modelica.SIunits.Area A_c = Modelica.Constants.pi * d_c * d_c / 8 "Area of semi-circular tube";    
  
  inner Modelica.SIunits.Area A_flow "Flow area of all channels";
  
  inner parameter Modelica.SIunits.Length length_cell = 1e-3 "length of a cell";
  
  inner Modelica.SIunits.Area A_stack "surface area of all cells in a stack";
    
  Modelica.SIunits.Area A_fc "Area of cold stream area";
  Modelica.SIunits.Area A_fh "Area of hot stream area";
  Modelica.SIunits.Area A_fmax "Area of maximum stream area comparing A_fc and A_fh: A_fmax = max(A_fc, A_fh)";  
  
  Modelica.SIunits.DynamicViscosity mu_c "average dynamic Viscosity in cold channel";
  Modelica.SIunits.DynamicViscosity mu_h "average dynamic Viscosity in hot channel"; 

  HXSegment segment[N_seg] "heat exchange segment";
  
  // length of one pipe in HeatExchanger
  Modelica.SIunits.Length length = length_cell * N_seg "length of a cell, unit m";
  
  inner Real fit_const_a "fitting constant a in Eq[3] of [kim, 2011] ";
  
  inner Real fit_const_b "fitting constant b in Eq[3] of [kim, 2011] ";
  
  inner Real fit_const_c "fitting constant c in Eq[4] of [kim, 2011] ";
  
  inner Real fit_const_d "fitting constant d in Eq[4] of [kim, 2011] ";
  
  //Modelica.SIunits.MassFlowRate G_c "Mass flux of each cold channel";  
  //Modelica.SIunits.MassFlowRate G_h "Mass flux of each cold channel";    
  
  Integer table_select;   
  
  model HXSegment
    "Segment in a HX, which contains a cold cell and hot cell"
    //extends BaseExchanger;
    //extends ExchangerContainer;
    
    replaceable package PBMedia = Steps.Media.SCO2;  
    
    replaceable Steps.Interfaces.PBFluidPort_a inlet_hot(redeclare package Medium = PBMedia) "Inlet port, previous component";
    replaceable Steps.Interfaces.PBFluidPort_b outlet_hot(redeclare package Medium = PBMedia) "Outlet port, next component";
    replaceable Steps.Interfaces.PBFluidPort_a inlet_cool(redeclare package Medium = PBMedia) "Recuperator inlet";
    replaceable Steps.Interfaces.PBFluidPort_b outlet_cool(redeclare package Medium = PBMedia) "Recuperator outlet";   
    
    outer Modelica.SIunits.Length t_wall "thickness of wall between two neighboring hot and cold";
    
    outer Modelica.SIunits.Area A_stack "surface area of all cells in a stack";
    
    outer Modelica.Blocks.Types.ExternalCombiTable1D table_th_inconel_750;
    
    outer String name_material;
    
    // Heat Flux
    inner Modelica.SIunits.HeatFlux q;  
    
    // wall thermal conductivity - determined by material of wall and local temperature
    Modelica.SIunits.ThermalConductivity k_wall;
    
    // overall Heat transfer coefficient
    Modelica.SIunits.CoefficientOfHeatTransfer U;      
  
    HXCell cell_cold(ByInlet = false);
    
    HXCell cell_hot(ByInlet = true);
    
    model HXCell    
      "One cell in a HX segment"
      
      replaceable package PBMedia = Steps.Media.SCO2; 
      
      replaceable Steps.Interfaces.PBFluidPort_a inlet(redeclare package Medium = Steps.Media.SCO2) "Inlet port, previous component";
      replaceable Steps.Interfaces.PBFluidPort_b outlet(redeclare package Medium = Steps.Media.SCO2) "Outlet port, next component";  
      
      outer Modelica.SIunits.Diameter d_h "Hydraulic Diameter";
      
      // Heat Flux
      //outer Modelica.SIunits.HeatFlux q;  
      
      outer Real fit_const_a "fitting constant a in Eq[3] of [kim, 2011] ";
      
      outer Real fit_const_b "fitting constant b in Eq[3] of [kim, 2011] ";
      
      outer Real fit_const_c "fitting constant c in Eq[4] of [kim, 2011] ";
      
      outer Real fit_const_d "fitting constant d in Eq[4] of [kim, 2011] ";
      
      outer Modelica.SIunits.Area A_flow "Flow area of all channels";
      
      outer Modelica.SIunits.Area A_c "Area of semi-circular tube"; 
      
      // length of this cell
      outer Modelica.SIunits.Length length_cell "unit m";  
      
      outer Integer N_channel "number of channels";     
      
      parameter Real id = 1 "id of the cell";   
            
      parameter Boolean debug_mode = false;
      
      parameter Boolean ByInlet = true "Flag for if the state is determined by Inlet, otherwise by outlet";
      
      // mass flow flux
      //Modelica.SIunits.MassFlowRate G;
      /*
      //Local temperature
      Modelica.SIunits.Temperature T;
      
      //Local pressure
      Modelica.SIunits.Pressure p;
      
      // Local velocity of fluid
      Modelica.SIunits.Velocity u;
      */
      // mass flow rate
      //Modelica.SIunits.MassFlowRate mdot;
      /*
      //local parameters of this cell listed as following
      //Local Dynamic Viscosity
      Modelica.SIunits.DynamicViscosity mu;
         
      // Local Conductivity
      Modelica.SIunits.ThermalConductivity k;
      
      // Local Reynolds Number
      Modelica.SIunits.ReynoldsNumber Re;
      
      // Local Density
      Modelica.SIunits.Density rho;
      
      // Local Nusselt Number
      Modelica.SIunits.NusseltNumber Nu;
      
      // Local PrandtlNumber
      //Modelica.SIunits.PrandtlNumber Pr;
      
      // local Thermal Conductance
      Modelica.SIunits.CoefficientOfHeatTransfer h;
      
      // Fanning Friction Factor - used to calculate pressure drop
      Real f;
      
      // Pressure drop
      Modelica.SIunits.PressureDifference dp;
      */
      // specific enthalpy to cal Heat flux
      Modelica.SIunits.SpecificEnthalpy h_mass;      
       
    algorithm   
      /*
      if ByInlet then
        p := inlet.p;
        T := inlet.T;
        mdot := inlet.m_flow;
      else
        p := outlet.p;
        T := outlet.T;  
        mdot := -outlet.m_flow;
      end if;
      
      G := mdot / N_channel / A_c;
      
      //Debug from this point
      mu := CP.PropsSI("V", "P", p, "T", T, PBMedia.mediumName); 
        
      k := CP.PropsSI("L", "P", p, "T", T, PBMedia.mediumName);  
      
      //h_mass := CP.PropsSI("H", "P", p, "T", T, PBMedia.mediumName);
      
      MyUtil.myAssert(debug = debug_mode, val_test = k, min = 0, max = 1e5, name_val = "k_c", val_ref = {id, T, p}, name_val_ref = {"id", "T", "P"});    
          
      Re := G * d_h / mu; 
        
      rho := CP.PropsSI("D", "P", p, "T", T, PBMedia.mediumName);
        
      u := mdot / A_flow / rho;    
      
      MyUtil.myAssert(debug = debug_mode, val_test = Re, min = 0, max = 1e6, name_val = "Re", name_val_ref = {"id", "G","d_h","mu"}, val_ref = {id, G, d_h, mu});
          
      Nu := 4.089 + fit_const_c * (Re ^ fit_const_d);
      
      h := 1.0; //Nu * k / d_h;
        
      f := (15.78 + fit_const_a * Re ^ fit_const_b ) / Re;       
      
      dp := 1.0; //f * length_cell * rho *  (u ^ 2) / d_h;    
      
      // no use of following parameters, use default value
      //Pr := 1.0;         
      */
      //mdot := 1.0;
      //q := 1.0;
      /*
      if ByInlet then        
        outlet.h_outflow := (inlet.h_outflow * 1.0 - 1.0) / 1.0;
        //outlet.T := CP.PropsSI("T", "P", outlet.p, "H", outlet.h_outflow , PBMedia.mediumName);
      else        
        inlet.h_outflow := (outlet.h_outflow * 1.0 + 1.0) / 1.0;
        //inlet.T := CP.PropsSI("T", "P", inlet.p, "H", inlet.h_outflow , PBMedia.mediumName);        
      end if; 
      */
      h_mass := 1.0;
    equation
      inlet.p - outlet.p = 0; //dp;
      inlet.T - outlet.T = 0; //dp;
      inlet.m_flow - outlet.m_flow = 0; //dp;
      //inlet.h_outflow - outlet.h_outflow = 0;
      inlet.h_outflow = inStream(outlet.h_outflow);
      
      if ByInlet then        
        outlet.h_outflow = (h_mass * 1.0 - 1.0) / 1.0;
        //outlet.T := CP.PropsSI("T", "P", outlet.p, "H", outlet.h_outflow , PBMedia.mediumName);
      else        
        inlet.h_outflow = (h_mass * 1.0 + 1.0) / 1.0;
        //inlet.T := CP.PropsSI("T", "P", inlet.p, "H", inlet.h_outflow , PBMedia.mediumName);        
      end if; 
      
      //(inlet.h_outflow - outlet.h_outflow) * 1.0 = 1.0; 
    		
    end HXCell;
   
  equation
  
    connect(inlet_cool, cell_cold.inlet);
    connect(cell_cold.outlet, outlet_cool);
    connect(inlet_hot, cell_hot.inlet);
    connect(cell_hot.outlet, outlet_hot);
   
    //connect(cell_hot.outlet_heat, cell_cold.inlet_heat);
      
  algorithm
  
    k_wall := 1.0;// MyUtil.thermal_conductivity(tableID = table_th_inconel_750, name = name_material, temperature = (cell_cold.T + cell_hot.T) / 2);
   
    U := 1.0; //1 / ( 1 / cell_hot.h + 1 / cell_cold.h + t_wall / k_wall);   
    /*
    if cell_hot.T > cell_cold.T then
      q := U * A_stack * (cell_hot.T - cell_cold.T);      
    else
      q := 0;
    end if;
    */
    
    q := 0.0;
  equation  
    
    inlet_cool.h_outflow = inStream(outlet_cool.h_outflow);
    inlet_hot.h_outflow = inStream(outlet_hot.h_outflow);
    
    //trivial equations, used to calculate mediums only. 
    /*   
    medium_cool_in.state = PBMedia.setState_pTX(p = inlet_cool.p, T = inlet_cool.T);
    medium_cool_out.state = PBMedia.setState_pTX(p = outlet_cool.p, T = outlet_cool.T);

    medium_hot_in.state = PBMedia.setState_pTX(p = inlet_hot.p, T = inlet_hot.T);
    medium_hot_out.state = PBMedia.setState_pTX(p = outlet_hot.p, T = outlet_hot.T);    
    
    inlet_cool_in.T = medium_cool_in.T;
    inlet_cool_in.p = medium_cool_in.T;
    inlet_cool_in.h_outflow = - medium_cool_in.h;
    inlet_cool_in.m_flow = - inlet_cool.m_flow;
    
    inlet_hot_in.T = medium_hot_in.T;
    inlet_hot_in.p = medium_hot_in.p;
    inlet_hot_in.h_outflow = - medium_hot_in.h;
    inlet_hot_in.m_flow = - inlet_hot.m_flow;
    
    outlet_cool_in.T = medium_cool_out.T;
    outlet_cool_in.p = medium_cool_out.p;
    outlet_cool_in.h_outflow = medium_cool_out.h;
    outlet_cool_in.m_flow = outlet_cool.m_flow;
    
    outlet_hot_in.T = medium_hot_out.T;
    outlet_hot_in.p = medium_hot_out.p;
    outlet_hot_in.h_outflow = medium_hot_out.h;
    outlet_hot_in.m_flow = outlet_hot.m_flow;
    */
  end HXSegment;
	
algorithm  

    when initial() then
      mu_c := CP.PropsSI("V", "P", p_cool, "T", T_cool_in, PBMedia.mediumName) "average dynamic Viscosity in cold channel";
      mu_h := CP.PropsSI("V", "P", p_hot, "T", T_hot_in, PBMedia.mediumName) "average dynamic Viscosity in hot channel"; 
      
      A_fc := m_dot_cool * d_h / mu_c / Re_design "Area of cold stream area";
      A_fh := m_dot_hot * d_h / mu_h /Re_design "Area of hot stream area";
      A_fmax := max(A_fc, A_fh) "Area of maximum stream area comparing A_fc and A_fh: A_fmax = max(A_fc, A_fh)";    
      
      t_wall := (2 - Modelica.Constants.pi  / 4) * (d_c / 2);
      
      N_channel := integer(A_fmax / A_c);
      
      A_flow := N_channel * A_c;
      
      A_stack := peri_c * length_cell * N_channel "surface area of all cells in a stack";    
      
      table_select := 4;
      // determine fitting constant by pitch and hydraulic diameter
      if (MyUtil.sameValue(pitch, 24.6 * 1e-3) and MyUtil.sameValue(d_h, 0.922 * 1e-3)) then
        table_select := 1;
        table_a := table_4a_a;
        table_b := table_4a_b;
        table_c := table_5a_c;
        table_d := table_5a_d;
      elseif (MyUtil.sameValue(pitch, 12.3 * 1e-3) and MyUtil.sameValue(d_h, 0.922 * 1e-3)) then
        table_select := 2;
        table_a := table_4b_a;
        table_b := table_4b_b;
        table_c := table_5b_c;
        table_d := table_5b_d;
      elseif (MyUtil.sameValue(pitch, 24.6 * 1e-3) and MyUtil.sameValue(d_h, 1.222 * 1e-3)) then
        table_select := 3;
        table_a := table_4c_a;
        table_b := table_4c_b;
        table_c := table_5c_c;
        table_d := table_5c_d;
      end if;
      
      fit_const_a := TB.CombiTable1D.getTableValue(table_a, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
      fit_const_b := TB.CombiTable1D.getTableValue(table_b, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
      fit_const_c := TB.CombiTable1D.getTableValue(table_c, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);
      fit_const_d := TB.CombiTable1D.getTableValue(table_d, icol = 1, u = Modelica.SIunits.Conversions.to_deg(phi), tableAvailable = 1.0);   
    
   end when;

equation  
  
  // connect all the segments within the heat exchanger, except for the end segment
  for i in 1 : N_seg loop
  
    // connect current segment's cool outlet with next segment's cool inlet
    if i <> N_seg then 
      connect(segment[i].outlet_cool, segment[i+1].inlet_cool);
    end if;
    
    // connect current segment's hot outlet with previous segment's hot inlet
    if i <> 1 then
      connect(segment[i].outlet_hot, segment[i-1].inlet_hot);
    end if;
        
  end for;
  
  // Now connect the end segement with my inlet and outlet
  connect(inlet_cool, segment[1].inlet_cool);
  connect(segment[N_seg].outlet_cool, outlet_cool);  
  
  connect(inlet_hot, segment[N_seg].inlet_hot);
  connect(segment[1].outlet_hot, outlet_hot);
/*
algorithm
  medium_hot_in.state := PBMedia.setState_pTX(p = inlet_hot.p, T = inlet_hot.T);
  medium_cool_out.state := PBMedia.setState_pTX(p = outlet_cool.p, T = outlet_cool.T);   
  
  medium_cool_in.state := PBMedia.setState_phX(p = segment[1].inlet_cool.p, h = segment[1].inlet_cool.h_outflow);       
  medium_hot_out.state := PBMedia.setState_phX(p = segment[N_seg].outlet_hot.p, h = segment[N_seg].outlet_hot.h_outflow);  
*/  
equation

  medium_hot_in.state = PBMedia.setState_pTX(p = inlet_hot.p, T = inlet_hot.T);
  medium_cool_out.state = PBMedia.setState_pTX(p = outlet_cool.p, T = outlet_cool.T);   
  
  medium_cool_in.state = PBMedia.setState_phX(p = segment[1].inlet_cool.p, h = segment[1].inlet_cool.h_outflow);       
  medium_hot_out.state = PBMedia.setState_phX(p = segment[N_seg].outlet_hot.p, h = segment[N_seg].outlet_hot.h_outflow);  

  //inlet_cool.T = medium_cool_in.T;
  //inlet_cool.p = medium_cool_in.p;
  //inlet_cool.h_outflow = segment[1].inlet.h_outflow;  
  //inlet_cool.h_outflow = inStream(outlet_cool.h_outflow);  
  //outlet_cool.h_outflow = - medium_cool_out.h;  
  inlet_cool.h_outflow = inStream(outlet_cool.h_outflow);
  
  /*
  inlet_cool_in.T = medium_cool_in.T;
  inlet_cool_in.p = medium_cool_in.T;
  inlet_cool_in.h_outflow = - medium_cool_in.h;
  inlet_cool_in.m_flow = - inlet_cool.m_flow;
  
  inlet_hot_in.T = medium_hot_in.T;
  inlet_hot_in.p = medium_hot_in.p;
  inlet_hot_in.h_outflow = - medium_hot_in.h;
  inlet_hot_in.m_flow = - inlet_hot.m_flow;
  
  outlet_cool_in.T = medium_cool_out.T;
  outlet_cool_in.p = medium_cool_out.p;
  outlet_cool_in.h_outflow = medium_cool_out.h;
  outlet_cool_in.m_flow = outlet_cool.m_flow;
  
  outlet_hot_in.T = medium_hot_out.T;
  outlet_hot_in.p = medium_hot_out.p;
  outlet_hot_in.h_outflow = medium_hot_out.h;
  outlet_hot_in.m_flow = outlet_hot.m_flow;
  */
  //outlet_cool.m_flow + inlet_cool.m_flow = 0;
  
  // for outlet_hot

  //outlet_hot.T = medium_hot_out.T;
  //outlet_hot.p = medium_hot_out.p;
  //outlet_hot.h_outflow = - segment[N_seg].outlet.h_outflow;    
  inlet_hot.h_outflow = inStream(outlet_hot.h_outflow);  
  //outlet_hot.m_flow + inlet_hot.m_flow = 0;
  
end MockPCHeatExchangerV2;
