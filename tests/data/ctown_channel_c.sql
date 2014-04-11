USE rtx_demo_db;
DROP TABLE IF EXISTS Msmts; 
DROP TABLE IF EXISTS Channels; 
# DROP TABLE IF EXISTS Distros;

CREATE TABLE Channels (  
			id		INTEGER AUTO_INCREMENT, 
			net_id_str	VARCHAR(255) NOT NULL,
			msmt_t	CHAR NOT NULL,
			l_lim	DOUBLE,
			r_lim	DOUBLE,
			info	VARCHAR(255),
			
			PRIMARY KEY (id)
			 
			); 
CREATE TABLE Msmts (
			id		INTEGER AUTO_INCREMENT, 
			time	DATETIME NOT NULL, 
			value	DOUBLE NOT NULL, 
			cid		INTEGER NOT NULL,
			PRIMARY KEY (id),
			CONSTRAINT fk_chnl FOREIGN KEY (cid) REFERENCES Channels(id) 
			);

# Test SCADA device suite for ctown complete network
# Control components
INSERT INTO Channels (net_id_str, msmt_t) VALUES 
	('PU1', 'C'), 	('PU2', 'C'), 	('V2', 'C'), 	('PU4', 'C'),
	('PU5', 'C'), 	('PU6', 'C'), 	('PU7', 'C'), 	('PU8', 'C'),
	('PU10', 'C'),  	('PU11', 'C')
	;
# Online sensors - tank levels and limits
INSERT INTO Channels (net_id_str, msmt_t, l_lim, r_lim) VALUES
	('R1', 'L', 0, 0), #reservior
	('T1', 'L', 71.5, 78), ('T2', 'L', 65, 70.9), ('T3', 'L', 112.9, 112.9+6.75), 
	('T4', 'L', 132.5, 132.5+4.7), ('T5', 'L', 105.8, 105.8+4.5), 
	('T6', 'L', 101.5, 101.5+5.5), ('T7', 'L', 102, 102+5) ;


# Online sensors - flowrates
INSERT INTO Channels (net_id_str, msmt_t, l_lim, r_lim) VALUES
	('PU1', 'Q', -200, 200), ('PU2', 'Q', -200, 200), ('PU3', 'Q', -200, 200), 
	('PU4', 'Q', -200, 200), ('PU5', 'Q', -200, 200), 
	('PU6', 'Q', -200, 200), ('PU7', 'Q', -200, 200), ('PU8', 'Q', -200, 200), 
	('PU9', 'Q', -200, 200), ('PU10', 'Q', -200, 200), ('PU11', 'Q', -200, 200),
	('P98', 'Q', -200, 200), ('P18', 'Q', -200, 200),('P946', 'Q', -200, 200),
	('P1016', 'Q', -200, 200), ('P138', 'Q', -200, 200),
	('P383', 'Q', -200, 200);

# Online sensors - pressures
INSERT INTO Channels (net_id_str, msmt_t) VALUES
	('J144', 'P'), ('J253', 'P'), ('J330', 'P'), ('J168', 'P'), ('J221', 'P'),
	('J28', 'P'); 


	