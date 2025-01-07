function [Filaments,Integrins,Ligands,FALconnections] = DeleteFALconnection(Field,Value,Filaments,Integrins,Ligands,FALconnections)


            switch Field
                    %...................................................................................
                    case 'FilamentName'  % Value = FilamentName
                            FALidx = find( FALconnections.FilamentName == Value ); 
                            if ~isempty(FALidx) % If this filament has a connection
                                 a = FALconnections.IntegrinIndex(FALidx,1);
                                 Integrins.AttachedFilamentName(a,1) = NaN;
                                 Integrins.AttachcedLigandIndex(a,1) = NaN;
                                 Integrins.ActiveStatus(a,1) = false;

                                 l = FALconnections.LigandIndex(FALidx,1);
                                 Ligands.AttachedFilamentName(l,1)  = NaN;
                                 Ligands.AttachedIntegrinIndex(l,1) = NaN;

                                 FALconnections.IntegrinIndex(FALidx,:) = [];
                                 FALconnections.FilamentName (FALidx,:) = [];
                                 FALconnections.MonomerIndex (FALidx,:) = [];
                                 FALconnections.LigandIndex  (FALidx,:) = [];
                            end


                    %...................................................................................
                    case 'FilamentNameAndMonomerIndex' % Needs to inputs for value: Value = [ FilamentName, MonomerIndex ] 
                            FALidx = find( FALconnections.FilamentName == Value(1) & ... % Find the integrin/ligand/filamement connection with this filament name
                                           FALconnections.MonomerIndex == Value(2) );    % and the integrin/ligand that is connected to the current monomer being depolymerized
                                % If so, reset integrin/ligand, and remove connection  
                                if ~isempty(FALidx) 
                                    a = FALconnections.IntegrinIndex(FALidx,1);
                                    Integrins.AttachedFilamentName (a,1) = NaN;
                                    Integrins.AttachedLigandIndex  (a,1) = NaN;
                                    Integrins.ActiveStatus         (a,1) = false;

                                    l = FALconnections.LigandIndex(FALidx);
                                    Ligands.AttachedFilamentName (l) = NaN;
                                    Ligands.AttachedIntegrinIndex(l) = NaN;

                                    FALconnections.IntegrinIndex(FALidx) = [];   
                                    FALconnections.LigandIndex  (FALidx) = [];
                                    FALconnections.FilamentName (FALidx) = []; 
                                    FALconnections.MonomerIndex (FALidx) = []; 
                                end


                    %...................................................................................
                    case 'IntegrinIndex'  % Value = IntegrinIndex
                                
                                for V = Value'
                                        FALidx = find( FALconnections.IntegrinIndex == V );
                                        if ~isempty(FALidx) 
                                            a = V;
                                            Integrins.AttachedFilamentName(a,1) = NaN;
                                            Integrins.AttachedLigandIndex (a,1) = NaN;
                                            Integrins.ActiveStatus        (a,1) = false;

                                            l = FALconnections.LigandIndex(FALidx);
                                            Ligands.AttachedFilamentName (l) = NaN;
                                            Ligands.AttachedIntegrinIndex(l) = NaN;

                                            FALconnections.IntegrinIndex(FALidx) = [];   
                                            FALconnections.LigandIndex  (FALidx) = [];
                                            FALconnections.FilamentName (FALidx) = []; 
                                            FALconnections.MonomerIndex (FALidx) = []; 
                                        end
                                end


                    %...................................................................................
                    case 'LigandIndex' % Value = LigandIndex
                            FALidx = find( FALconnections.LigandIndex == Value );
                            if ~isempty(FALidx) 
                                a = FALconnections.IntegrinIndex(FALidx,1);
                                Integrins.AttachedFilamentName(a,1) = NaN;
                                Integrins.AttachedLigandIndex (a,1) = NaN;
                                Integrins.ActiveStatus        (a,1) = false;

                                l = FALconnections.LigandIndex(FALidx);
                                Ligands.AttachedFilamentName (l) = NaN;
                                Ligands.AttachedIntegrinIndex(l) = NaN;

                                FALconnections.IntegrinIndex(FALidx) = [];   
                                FALconnections.LigandIndex  (FALidx) = [];
                                FALconnections.FilamentName (FALidx) = []; 
                                FALconnections.MonomerIndex (FALidx) = []; 
                            end
                            
                            
                    %...................................................................................
            end



end
