within Steps.Components;

model PCHeatExchangerV2
  "Printed Circuit based Heat Exchanger - version 2"
  extends BaseExchanger;  

  import CP = Steps.Utilities.CoolProp;  
  import UTIL = Modelica.Utilities;
  import MyUtil = Steps.Utilities.Util;
  
  replaceable package PBMedia = Steps.Media.SCO2;   
  
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
  
  parameter Boolean debug_mode = false;  
  
  
protected  
  
  inner KimCorrelations kim_cor(phi = phi, pitch = pitch, d_h = d_h);
  
  inner parameter String name_material = "inconel 750";  
  
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
  
  //Modelica.SIunits.MassFlowRate G_c "Mass flux of each cold channel";  
  //Modelica.SIunits.MassFlowRate G_h "Mass flux of each cold channel";    
  
  Integer table_select;   
  
  model HXSegment
    "Segment in a HX, which contains a cold cell and hot cell"
    extends BaseExchanger;
    replaceable package PBMedia = Steps.Media.SCO2;  
   
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
      
      parameter Real id = 1 "id of the cell";     
      
      outer Modelica.SIunits.Diameter d_h "Hydraulic Diameter";
      
      // Heat Flux
      outer Modelica.SIunits.HeatFlux q;  
      
      outer KimCorrelations kim_cor;      
     
      outer Modelica.SIunits.Area A_flow "Flow area of all channels";
      
      outer Modelica.SIunits.Area A_c "Area of semi-circular tube"; 
      
      // length of this cell
      outer Modelica.SIunits.Length length_cell "unit m";  
      
      outer Integer N_channel "number of channels";
      
      Modelica.SIunits.MassFlowRate G;
      
      parameter Boolean debug_mode = false;
      
      parameter Boolean ByInlet = true "Flag for if the state is determined by Inlet, otherwise by outlet";
      
      //Local temperature
      Modelica.SIunits.Temperature T;
      
      //Local pressure
      Modelica.SIunits.Pressure p;
      
      // Local velocity of fluid
      Modelica.SIunits.Velocity u;
      
      // mass flow rate
      Modelica.SIunits.MassFlowRate mdot;
      
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
      
      // specific enthalpy to cal Heat flux
      //Modelica.SIunits.SpecificEnthalpy h_mass;      
       
    algorithm   
    
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
          
      Nu := 4.089 + kim_cor.c * (Re ^ kim_cor.d);
       
      h := Nu * k / d_h;
        
      f := (15.78 + kim_cor.a * Re ^ kim_cor.b ) / Re;       
      
      dp := f * length_cell * rho *  (u ^ 2) / d_h;    
      
      // no use of following parameters, use default value
      //Pr := 1.0;         
    
      if ByInlet then        
        //outlet.h_outflow := (inlet.h_outflow * mdot - q) / mdot;
        outlet.T := CP.PropsSI("T", "P", outlet.p, "H", outlet.h_outflow , PBMedia.mediumName);
      else        
        //inlet.h_outflow := (outlet.h_outflow * mdot + q) / mdot;
        inlet.T := CP.PropsSI("T", "P", inlet.p, "H", inlet.h_outflow , PBMedia.mediumName);        
      end if; 
    
    equation
      inlet.p - outlet.p = dp;
    
      //inlet.h_outflow = inStream(outlet.h_outflow);
      
      (inlet.h_outflow - outlet.h_outflow) * mdot = q; 
    		
    end HXCell;
    /*
    model HXCellCold
     extends HXCell;
    
     Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a inlet_heat "Inlet of heat flow";
    
    
    end HXCellCold;
    
    model HXCellHot
      extends HXCell;
    
      Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b outlet_heat "Outlet of heat flow";
    
    end HXCellHot; 		
    */  
  
  equation
  
    connect(inlet_cool, cell_cold.inlet);
    connect(cell_cold.outlet, outlet_cool);
    connect(inlet_hot, cell_hot.inlet);
    connect(cell_hot.outlet, outlet_hot);
   
    //connect(cell_hot.outlet_heat, cell_cold.inlet_heat);
      
  algorithm
  
    k_wall := MyUtil.thermal_conductivity(tableID = table_th_inconel_750, name = name_material, temperature = (cell_cold.T + cell_hot.T) / 2);
   
    U := 1 / ( 1 / cell_hot.h + 1 / cell_cold.h + t_wall / k_wall);   
   
    if cell_hot.T > cell_cold.T then
      q := U * A_stack * (cell_hot.T - cell_cold.T);      
    else
      q := 0;
    end if;
    
  equation  
    
    //trivial equations, used to calculate mediums only. 
        
    medium_cool_in.state = PBMedia.setState_pTX(p = inlet_cool.p, T = inlet_cool.T);
    medium_cool_out.state = PBMedia.setState_pTX(p = outlet_cool.p, T = outlet_cool.T);

    medium_hot_in.state = PBMedia.setState_pTX(p = inlet_hot.p, T = inlet_hot.T);
    medium_hot_out.state = PBMedia.setState_pTX(p = outlet_hot.p, T = outlet_hot.T);    
    
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

algorithm
  medium_hot_in.state := PBMedia.setState_pTX(p = inlet_hot.p, T = inlet_hot.T);
  medium_cool_out.state := PBMedia.setState_pTX(p = outlet_cool.p, T = outlet_cool.T);   
  
  medium_cool_in.state := PBMedia.setState_phX(p = segment[1].inlet_cool.p, h = segment[1].inlet_cool.h_outflow);       
  medium_hot_out.state := PBMedia.setState_phX(p = segment[N_seg].outlet_hot.p, h = segment[N_seg].outlet_hot.h_outflow);  
  
equation
  //inlet_cool.T = medium_cool_in.T;
  //inlet_cool.p = medium_cool_in.p;
  //inlet_cool.h_outflow = segment[1].inlet.h_outflow;  
  //inlet_cool.h_outflow = inStream(outlet_cool.h_outflow);  
  //outlet_cool.h_outflow = - medium_cool_out.h;  
  inlet_cool.h_outflow = - outlet_cool.h_outflow;
  //outlet_cool.m_flow + inlet_cool.m_flow = 0;
  
  // for outlet_hot

  //outlet_hot.T = medium_hot_out.T;
  //outlet_hot.p = medium_hot_out.p;
  //outlet_hot.h_outflow = - segment[N_seg].outlet.h_outflow;    
  inlet_hot.h_outflow = inStream(outlet_hot.h_outflow);  
  //outlet_hot.m_flow + inlet_hot.m_flow = 0;
  
end PCHeatExchangerV2;
