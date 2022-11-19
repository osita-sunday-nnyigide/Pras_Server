/*******************************************************************************************************************************
This file is a part of the Protein Repair and Analysis Server written by Osita S. Nnyigide

See the included LICENSE file
********************************************************************************************************************************/

#include "SecondaryStructure.hpp"

  void AnalyzeProtein::SecondaryStructure(char ftype)
  {
  	cout<<"**Now assigning secondary structure elements...\n"<<endl;
  	VS  cid, SSA; VVS resn, resno, dtype, tSSA; VI cnt;
  	VVF ahlx, amdih, _3hlx, btn, phid, psid, phlx;
    if (ftype == 'b')      getChains("out.pdb", " ", " ", " ", " ");
    else if (ftype == 'f') getMMCIFChains("out.cif", " ", " ", " ");
    resn  = resName();     
    resno = resNumber();     
    cid   = chainID();
    dtype = dihedralType();
    ahlx  = alpha_helix();   
    amdih = amideDihedral();
    _3hlx = _310_helix(); 
    btn    = beta_turn();      
    phid = phiDihedral();
    psid  = psiDihedral(); 
    phlx  = pi_helix();
    tSSA = {};
		for (int k= 0; k<cid.size(); k++)
		{
			SSA = {};
			for (int l= 0; l<ahlx[k].size(); l++)
			{
				if    ((amdih[k][l] > -10.0 && amdih[k][l] < 45.0 && phlx[k][l-1] < 3.5 && amdih[k][(l+1)%amdih[k].size()] > -10.
					  &&amdih[k][(l+1)%amdih[k].size()]<45.&&amdih[k][(l+2)%amdih[k].size()]>-10.&&amdih[k][(l+2)%amdih[k].size()]<45.)
					  &&((phlx[k][l-1] < ahlx[k][l-1] + 0.1)||(phlx[k][l] < ahlx[k][l] + 0.1 && phlx[k][(l+1)%phlx[k].size()] <
					  ahlx[k][(l+1)%ahlx[k].size()] + 0.1))) 
					{
					  SSA.push_back("Pi_helix");
					}
				else if  ((amdih[k][l] > -10.0 && amdih[k][l] < 45.0 && ahlx[k][l-1] < 3.5 && amdih[k][(l+1)%amdih[k].size()]>-10.
				         && amdih[k][(l+1)%amdih[k].size()]<45.0&&amdih[k][(l+2)%amdih[k].size()]>-10.0&&amdih[k][(l+2)%amdih[k].size()]<45.0)
			             && ((ahlx[k][l-1] < _3hlx[k][l-1] + 0.1) || (ahlx[k][l] < _3hlx[k][l] + 0.1 &&
				         ahlx[k][(l+1)%ahlx[k].size()] < _3hlx[k][(l+1)%_3hlx[k].size()] + 0.1))) 
					{
						 SSA.push_back("Alpha_helix");
					}
				else if  ((amdih[k][l] > -10. && amdih[k][l] < 45. && _3hlx[k][l-1] < 3.5 && amdih[k][(l+1)%amdih[k].size()]>-10.
				         &&amdih[k][(l+1)%amdih[k].size()]<45.&&amdih[k][(l+2)%amdih[k].size()]>-10.&&amdih[k][(l+2)%amdih[k].size()]<45)
			             && ((_3hlx[k][l-1] < ahlx[k][l-1] + 0.1) || (_3hlx[k][l] < ahlx[k][l] + 0.1 &&
				         _3hlx[k][(l+1)%_3hlx[k].size()] < ahlx[k][(l+1)%ahlx[k].size()] + 0.1))) 
					{
						 SSA.push_back("310_helix");
					} 
				else
					{
						SSA.push_back("Coil");
					}
			}
			tSSA.push_back(SSA);
		}
    for (int i= 0; i<tSSA.size(); i++)
	    {
	    	for (int j= 0; j<tSSA[i].size(); j++)
		    	if ((tSSA[i][j] == "Pi_helix") && (tSSA[i][(j+1)%tSSA[i].size()] != "Pi_helix") && (VS {tSSA[i].begin()+j, tSSA[i].end()}).size()>4)
			    	{
			    		cnt.push_back(j);
			    	}
			  for (auto n:cnt)
			  		{
			  			tSSA[i][n+1] = "Pi_helix";
			    		tSSA[i][n+2] = "Pi_helix";
			    		tSSA[i][n+3] = "Pi_helix";
			    		tSSA[i][n+4] = "Pi_helix";
			    		cnt.clear();
			  		} 
	    	for (int j= 0; j<tSSA[i].size(); j++)
		    	 if ((tSSA[i][j] == "Alpha_helix")&&(tSSA[i][(j+1)%tSSA[i].size()]!="Pi_helix" && tSSA[i][(j+1)%tSSA[i].size()]!="Alpha_helix")&&
		    		    (VS {tSSA[i].begin()+j, tSSA[i].end()}).size()>3)
			    	{
			    		cnt.push_back(j);
			    	}
			  for (auto n:cnt)
			  		{
			  			tSSA[i][n+1] = "Alpha_helix";
			    		tSSA[i][n+2] = "Alpha_helix";
			    		tSSA[i][n+3] = "Alpha_helix";
			    		cnt.clear();
			  		} 
	    for (int j= 0; j<tSSA[i].size(); j++)
		    	 if ((tSSA[i][j] == "310_helix") && (tSSA[i][(j+1)%tSSA[i].size()] != "310_helix") && (VS {tSSA[i].begin()+j, tSSA[i].end()}).size()>2)
			    	{
							cnt.push_back(j);
			    	}
			  for (auto n:cnt)
			  		{
			  			tSSA[i][n+1] = "310_helix";
			    		tSSA[i][n+2] = "310_helix";
			    		cnt.clear();
			  		} 
			for (int j= 0; j<ahlx[i].size(); j++)
				  if  ((((amdih[i][j]>173&&amdih[i][j]<181)||(amdih[i][j]<-140))&&((amdih[i][(j+1)%amdih[i].size()]>=173&&amdih[i][(j+1)%amdih[i].size()]<181) 
					   ||(amdih[i][(j+1)%amdih[i].size()]<=-140)))||(phid[i][j]<-120&&psid[i][j]>120&&tSSA[i][j]!="Alpha_helix" && tSSA[i][j] != "310_helix" &&
					   tSSA[i][j]!="Pi_helix")&&(tSSA[i][(j+1)%tSSA[i].size()]!="Alpha_helix"&&tSSA[i][(j+1)%tSSA[i].size()]!="Pi_helix"&&tSSA[i][(j+1)%tSSA[i].size()]
					   !="310_helix") && (VS {tSSA[i].begin()+j, tSSA[i].end()}).size()>2)
						{
							cnt.push_back(j);				
						}
			  for (auto n:cnt)
			  		{
					    tSSA[i][n]   = "Strand";
					    tSSA[i][n+1] = "Strand";	
			    		cnt.clear();
			  		}
			for (int j= 0; j<ahlx[i].size(); j++)
				  if ((_3hlx[i][j]<=3.5)&&(tSSA[i][j]!="Alpha_helix"&&tSSA[i][j]!="310_helix"&&tSSA[i][j]!="Pi_helix")&& (tSSA[i][(j+1)%tSSA[i].size()]
				    		!="Alpha_helix" && tSSA[i][(j+1)%tSSA[i].size()]!="Pi_helix"&& tSSA[i][(j+1)%tSSA[i].size()] != "310_helix") &&  (VS {tSSA[i].begin()+j, 
				    		tSSA[i].end()}).size()>2)
						{
							cnt.push_back(j);				
						}
			  for (auto n:cnt)
			  		{
					    tSSA[i][n]   = "Turn_hb";
					    tSSA[i][n+1] = "Turn_hb";	
			    		cnt.clear();
			  		}
			for (int j= 0; j<ahlx[i].size(); j++)
					 if  ((btn[i][j] <= 4.7) && (tSSA[i][j]!= "Alpha_helix" && tSSA[i][j] != "310_helix" && tSSA[i][j] != "Pi_helix" && tSSA[i][j]
							 != "Turn_hb" && tSSA[i][j] != "Strand") && (tSSA[i][(j+1)%tSSA[i].size()] != "Alpha_helix" && tSSA[i][(j+1)%tSSA[i].size()]
							 != "Pi_helix" && tSSA[i][(j+1)%tSSA[i].size()]!= "310_helix") &&  (VS {tSSA[i].begin()+j, tSSA[i].end()}).size()>2)
						{
							cnt.push_back(j);				
						}
			  for (auto n:cnt)
			  		{
					    tSSA[i][n]   = "Turn_nhb";
					    tSSA[i][n+1] = "Turn_nhb";	
			    		cnt.clear();
			  		}
			for (int j= 0; j<ahlx[i].size(); j++)
					 if ((amdih[i][j] < -150 && amdih[i][j] < -80) && (phid[i][j] < -65 && phid[i][j] > -85 && psid[i][j] < 175 && psid[i][j] > 140)
						    && (phid[i][j-1] < -65 && phid[i][j-1] > -85 && psid[i][j] < 175 && psid[i][j] > 140) && (tSSA[i][j] != "Alpha_helix")
						    && (tSSA[i][j] != "310_helix") && (tSSA[i][j] != "Pi_helix") && (tSSA[i][j] != "Turn_hb") && (tSSA[i][j] != "Turn_nhb")
						    && (tSSA[i][j] != "Strand") &&  (VS {tSSA[i].begin()+j, tSSA[i].end()}).size()>1)
						{
						    tSSA[i][j] = "Polyp_helix";
						}
			float total_coil       =  count(tSSA[i].begin(), tSSA[i].end(), "Coil");
			float total_nhb_turn   =  count(tSSA[i].begin(), tSSA[i].end(), "Turn_nhb");
			float total_hb_turn    =  count(tSSA[i].begin(), tSSA[i].end(), "Turn_hb");
			float total_strand     =  count(tSSA[i].begin(), tSSA[i].end(), "Strand");
			float total_alphahelix =  count(tSSA[i].begin(), tSSA[i].end(), "Alpha_helix");
			float total_310helix   =  count(tSSA[i].begin(), tSSA[i].end(), "310_helix");
			float total_pihelix    =  count(tSSA[i].begin(), tSSA[i].end(), "Pi_helix");
			float total_pphelix    =  count(tSSA[i].begin(), tSSA[i].end(), "Polyp_helix");
			float pcent_ahelix     =  (total_alphahelix/tSSA[i].size())*100.;
			float pcent_phelix     =  (total_pihelix/tSSA[i].size())*100.;
			float pcent_310helix   =  (total_310helix/tSSA[i].size())*100.;
			float pcent_pphelix    =  (total_pphelix/tSSA[i].size())*100.;
			float pcent_strand     =  (total_strand/tSSA[i].size())*100.;	
			float pcent_turn       =  ((total_nhb_turn+total_hb_turn)/tSSA[i].size())*100.;	
		  float pcent_coil       =  (total_coil/tSSA[i].size())*100.;									
			string pdb_code        =  "Chain";
			ofstream ssa;
		    ssa.open("sec_strc_table_chain"+to_string(i+1)+".txt");
		    ssa<<setiosflags(ios::fixed)<<setprecision(3);
		    ssa<<"When using files generated by this program in a publication, please cite this program as O.S. Nnyigide, T.O. Nnyigide"
		               " S.G. Lee, K. Hyun. Protein Repair and Analysis Server: A Web\nServer to Repair PDB Structures, Add Missing Heavy"
		    " Atoms and Hydrogen Atoms, and Assign Secondary Structures by Amide Interactions. J. Chem. Inf. Model., 2022, 62, 4232–4246\n";
		    ssa<<string(183,'*')<<"\n";ssa<<"Secondary Structure Assignment by Amide-Amide Interactions of the Backbone"
		            " This is the summary of the measurements used to assign the secondary structure elements\n";ssa<<string(183,'*')<<endl;

			ssa<<"  Residues in PDB    Amide_dihedral   Phi_dihedral   Psi_dihedral     O(i)->N(i+4)    O(i)->N(i+3)"
					"    O(i)->N(i+5)    C(i)->N(i+3)   Dihedral_type   Secondary structure"<<endl;

			for (int k= 0; k<tSSA[i].size(); k++)
			{
               ssa<<setw(15)<<pdb_code+cid[i]+":"+resn[i][k]+":"+resno[i][k]   <<setw(16)<<amdih[i][k]   <<setw(16)<<phid[i][k]   <<setw(16)<<psid[i][k]
               <<setw(16)<<ahlx[i][k]<<setw(16)<<_3hlx[i][k]<<setw(16)<<phlx[i][k]<<setw(16)<<btn[i][k]<<setw(16)<<dtype[i][k]<<setw(16)<<tSSA[i][k]<<endl;
			}
			ofstream raman;
		    raman.open("rama_ssa_chain"+to_string(i+1)+".txt");
		    for (int k= 0; k<tSSA[i].size(); k++)
		    {
		    	raman<<dtype[i][k]<<" "<<phid[i][k]<<" "<<psid[i][k]<<endl;
		    }
		    if (pcent_ahelix > 0.001)raman<<"Alpha_helix"<<" "<<pcent_ahelix<<endl;
		    if (pcent_phelix > 0.001)raman<<"Pi_helix"<<" "<<pcent_phelix<<endl;
		    if (pcent_310helix > 0.001)raman<<"310_helix"<<" "<<pcent_310helix<<endl;
		    if (pcent_pphelix > 0.001)raman<<"Polyp_helix"<<" "<<pcent_pphelix<<endl;
		    if (pcent_strand > 0.001)raman<<"Strand"<<" "<<pcent_strand<<endl;
		    if (pcent_turn > 0.001)raman<<"Turn"<<" "<<pcent_turn<<endl;
		    if (pcent_coil > 0.001)raman<<"Coil"<<" "<<pcent_coil<<endl;
	    }
  }
