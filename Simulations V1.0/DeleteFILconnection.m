function [Filaments,Integrins,Ligands,FILconnections] = DeleteFILconnection(Field,Value,Filaments,Integrins,Ligands,FILconnections)


            switch Field
                    %...................................................................................
                    case 'FilamentName'  % Value = FilamentName
                            FILidx = find( FILconnections.FilamentName == Value ); 
                            if ~isempty(FILidx) % If this filament has a connection
                                 a = FILconnections.IntegrinIndex(FILidx,1);
                                 Integrins.AttachedFilamentName(a,1) = NaN;
                                 Integrins.AttachcedLigandIndex(a,1) = NaN;
                                 Integrins.ActiveStatus(a,1) = false;

                                 l = FILconnections.LigandIndex(FILidx,1);
                                 Ligands.AttachedFilamentName(l,1)  = NaN;
                                 Ligands.AttachedIntegrinIndex(l,1) = NaN;

                                 FILconnections.IntegrinIndex(FILidx,:) = [];
                                 FILconnections.FilamentName (FILidx,:) = [];
                                 FILconnections.MonomerIndex (FILidx,:) = [];
                                 FILconnections.LigandIndex  (FILidx,:) = [];
                            end


                    %...................................................................................
                    case 'FilamentNameAndMonomerIndex' % Needs to inputs for value: Value = [ FilamentName, MonomerIndex ] 
                            FILidx = find( FILconnections.FilamentName == Value(1) & ... % Find the integrin/ligand/filamement connection with this filament name
                                           FILconnections.MonomerIndex == Value(2) );    % and the integrin/ligand that is connected to the current monomer being depolymerized
                                % If so, reset integrin/ligand, and remove connection  
                                if ~isempty(FILidx) 
                                    a = FILconnections.IntegrinIndex(FILidx,1);
                                    Integrins.AttachedFilamentName (a,1) = NaN;
                                    Integrins.AttachedLigandIndex  (a,1) = NaN;
                                    Integrins.ActiveStatus         (a,1) = false;

                                    l = FALconnections.LigandIndex(FILidx);
                                    Ligands.AttachedFilamentName (l) = NaN;
                                    Ligands.AttachedIntegrinIndex(l) = NaN;

                                    FILconnections.IntegrinIndex(FILidx) = [];   
                                    FILconnections.LigandIndex  (FILidx) = [];
                                    FILconnections.FilamentName (FILidx) = []; 
                                    FILconnections.MonomerIndex (FILidx) = []; 
                                end


                    %...................................................................................
                    case 'IntegrinIndex'  % Value = IntegrinIndex
                                
                                for V = Value'
                                        FILidx = find( FILconnections.IntegrinIndex == V );
                                        if ~isempty(FILidx) 
                                            a = V;
                                            Integrins.AttachedFilamentName(a,1) = NaN;
                                            Integrins.AttachedLigandIndex (a,1) = NaN;
                                            Integrins.ActiveStatus        (a,1) = false;

                                            l = FALconnections.LigandIndex(FILidx);
                                            Ligands.AttachedFilamentName (l) = NaN;
                                            Ligands.AttachedIntegrinIndex(l) = NaN;

                                            FILconnections.IntegrinIndex(FILidx) = [];   
                                            FILconnections.LigandIndex  (FILidx) = [];
                                            FILconnections.FilamentName (FILidx) = []; 
                                            FILconnections.MonomerIndex (FILidx) = []; 
                                        end
                                end


                    %...................................................................................
                    case 'LigandIndex' % Value = LigandIndex
                            FILidx = find( FILconnections.LigandIndex == Value );
                            if ~isempty(FILidx) 
                                a = FILconnections.IntegrinIndex(FALidx,1);
                                Integrins.AttachedFilamentName(a,1) = NaN;
                                Integrins.AttachedLigandIndex (a,1) = NaN;
                                Integrins.ActiveStatus        (a,1) = false;

                                l = FILconnections.LigandIndex(FALidx);
                                Ligands.AttachedFilamentName (l) = NaN;
                                Ligands.AttachedIntegrinIndex(l) = NaN;

                                FILconnections.IntegrinIndex(FALidx) = [];   
                                FILconnections.LigandIndex  (FALidx) = [];
                                FILconnections.FilamentName (FALidx) = []; 
                                FILconnections.MonomerIndex (FALidx) = []; 
                            end
                            
                            
                    %...................................................................................
            end



end
