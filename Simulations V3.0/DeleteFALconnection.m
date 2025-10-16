function [Filaments,Adhesions,Ligands,FALconnections] = DeleteFALconnection(Field,Value,Filaments,Adhesions,Ligands,FALconnections)


            switch Field
                    %...................................................................................
                    case 'FilamentName'  % Value = FilamentName
                            FALidx = find( FALconnections.FilamentName == Value ); 
                            if ~isempty(FALidx) % If this filament has a connection
                                 a = FALconnections.AdhesionIndex(FALidx,1);
                                 Adhesions.AttachedFilamentName(a,1) = NaN;
                                 Adhesions.AttachcedLigandIndex(a,1) = NaN;
                                 Adhesions.ActiveStatus(a,1) = false;

                                 l = FALconnections.LigandIndex(FALidx,1);
                                 Ligands.AttachedFilamentName(l,1)  = NaN;
                                 Ligands.AttachedAdhesionIndex(l,1) = NaN;

                                 FALconnections.AdhesionIndex(FALidx,:) = [];
                                 FALconnections.FilamentName (FALidx,:) = [];
                                 FALconnections.MonomerIndex (FALidx,:) = [];
                                 FALconnections.LigandIndex  (FALidx,:) = [];
                            end


                    %...................................................................................
                    case 'FilamentNameAndMonomerIndex' % Needs to inputs for value: Value = [ FilamentName, MonomerIndex ] 
                            FALidx = find( FALconnections.FilamentName == Value(1) & ... % Find the adhesion/ligand/filamement connection with this filament name
                                           FALconnections.MonomerIndex == Value(2) );    % and the adhesion/ligand that is connected to the current monomer being depolymerized
                                % If so, reset adhesion/ligand, and remove connection  
                                if ~isempty(FALidx) 
                                    a = FALconnections.AdhesionIndex(FALidx,1);
                                    Adhesions.AttachedFilamentName (a,1) = NaN;
                                    Adhesions.AttachedLigandIndex  (a,1) = NaN;
                                    Adhesions.ActiveStatus         (a,1) = false;

                                    l = FALconnections.LigandIndex(FALidx);
                                    Ligands.AttachedFilamentName (l) = NaN;
                                    Ligands.AttachedAdhesionIndex(l) = NaN;

                                    FALconnections.AdhesionIndex(FALidx) = [];   
                                    FALconnections.LigandIndex  (FALidx) = [];
                                    FALconnections.FilamentName (FALidx) = []; 
                                    FALconnections.MonomerIndex (FALidx) = []; 
                                end


                    %...................................................................................
                    case 'AdhesionIndex'  % Value = AdhesionIndex
                                
                                for V = Value'
                                        FALidx = find( FALconnections.AdhesionIndex == V );
                                        if ~isempty(FALidx) 
                                            a = V;
                                            Adhesions.AttachedFilamentName(a,1) = NaN;
                                            Adhesions.AttachedLigandIndex (a,1) = NaN;
                                            Adhesions.ActiveStatus        (a,1) = false;

                                            l = FALconnections.LigandIndex(FALidx);
                                            Ligands.AttachedFilamentName (l) = NaN;
                                            Ligands.AttachedAdhesionIndex(l) = NaN;

                                            FALconnections.AdhesionIndex(FALidx) = [];   
                                            FALconnections.LigandIndex  (FALidx) = [];
                                            FALconnections.FilamentName (FALidx) = []; 
                                            FALconnections.MonomerIndex (FALidx) = []; 
                                        end
                                end


                    %...................................................................................
                    case 'LigandIndex' % Value = LigandIndex
                            FALidx = find( FALconnections.LigandIndex == Value );
                            if ~isempty(FALidx) 
                                a = FALconnections.AdhesionIndex(FALidx,1);
                                Adhesions.AttachedFilamentName(a,1) = NaN;
                                Adhesions.AttachedLigandIndex (a,1) = NaN;
                                Adhesions.ActiveStatus        (a,1) = false;

                                l = FALconnections.LigandIndex(FALidx);
                                Ligands.AttachedFilamentName (l) = NaN;
                                Ligands.AttachedAdhesionIndex(l) = NaN;

                                FALconnections.AdhesionIndex(FALidx) = [];   
                                FALconnections.LigandIndex  (FALidx) = [];
                                FALconnections.FilamentName (FALidx) = []; 
                                FALconnections.MonomerIndex (FALidx) = []; 
                            end
                            
                            
                    %...................................................................................
            end



end