function [Filaments,Integrins,Ligands,FILconnections] = DeleteFILconnection(Field,Value,Filaments,Integrins,Ligands,FILconnections)


            switch Field
                    %...................................................................................
                    case 'FilamentName'  % Value = FilamentName
                            FILidx = find( FILconnections.FilamentName == Value ); 
                            if ~isempty(FALidx) % If this filament has a connection
                                 a = FILconnections.IntegrinIndex(FALidx,1);
                                 Integrins.AttachedFilamentName(a,1) = NaN;
                                 Integrins.AttachcedLigandIndex(a,1) = NaN;
                                 Integrins.ActiveStatus(a,1) = false;

                                 l = FILconnections.LigandIndex(FALidx,1);
                                 Ligands.AttachedFilamentName(l,1)  = NaN;
                                 Ligands.AttachedIntegrinIndex(l,1) = NaN;

                                 FILconnections.IntegrinIndex(FALidx,:) = [];
                                 FILconnections.FilamentName (FALidx,:) = [];
                                 FILconnections.MonomerIndex (FALidx,:) = [];
                                 FILconnections.LigandIndex  (FALidx,:) = [];
                            end


                    %...................................................................................
                    case 'FilamentNameAndMonomerIndex' % Needs to inputs for value: Value = [ FilamentName, MonomerIndex ] 
                            FILidx = find( FILconnections.FilamentName == Value(1) & ... % Find the integrin/ligand/filamement connection with this filament name
                                           FILconnections.MonomerIndex == Value(2) );    % and the integrin/ligand that is connected to the current monomer being depolymerized
                                % If so, reset integrin/ligand, and remove connection  
                                if ~isempty(FILidx) 
                                    a = FILconnections.IntegrinIndex(FALidx,1);
                                    Integrins.AttachedFilamentName (a,1) = NaN;
                                    Integrins.AttachedLigandIndex  (a,1) = NaN;
                                    Integrins.ActiveStatus         (a,1) = false;

                                    l = FALconnections.LigandIndex(FALidx);
                                    Ligands.AttachedFilamentName (l) = NaN;
                                    Ligands.AttachedIntegrinIndex(l) = NaN;

                                    FILconnections.IntegrinIndex(FALidx) = [];   
                                    FILconnections.LigandIndex  (FALidx) = [];
                                    FILconnections.FilamentName (FALidx) = []; 
                                    FILconnections.MonomerIndex (FALidx) = []; 
                                end


                    %...................................................................................
                    case 'IntegrinIndex'  % Value = IntegrinIndex
                                
                                for V = Value'
                                        FALidx = find( FILconnections.IntegrinIndex == V );
                                        if ~isempty(FILidx) 
                                            a = V;
                                            Integrins.AttachedFilamentName(a,1) = NaN;
                                            Integrins.AttachedLigandIndex (a,1) = NaN;
                                            Integrins.ActiveStatus        (a,1) = false;

                                            l = FALconnections.LigandIndex(FALidx);
                                            Ligands.AttachedFilamentName (l) = NaN;
                                            Ligands.AttachedIntegrinIndex(l) = NaN;

                                            FILconnections.IntegrinIndex(FALidx) = [];   
                                            FILconnections.LigandIndex  (FALidx) = [];
                                            FILconnections.FilamentName (FALidx) = []; 
                                            FILconnections.MonomerIndex (FALidx) = []; 
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
