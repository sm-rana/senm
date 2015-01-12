set search_path to scada,public;
DROP TABLE IF EXISTS msmts; 
DROP TABLE IF EXISTS channels; 
-- DROP TABLE IF EXISTS Distros;

CREATE TABLE channels (  
			id		SERIAL, 
			net_id_str	VARCHAR(255) NOT NULL,
			msmt_t	CHAR NOT NULL,
			l_lim	DOUBLE PRECISION,
			r_lim	DOUBLE PRECISION,
			info	VARCHAR(255),
			
			PRIMARY KEY (id)
			 
			); 
create index on channels (id);

CREATE TABLE msmts (
			id		SERIAL, 
			time	TIMESTAMP NOT NULL, 
			value	DOUBLE PRECISION NOT NULL, 
			cid		INTEGER NOT NULL,
			PRIMARY KEY (id),
			CONSTRAINT fk_chnl FOREIGN KEY (cid) REFERENCES Channels(id) 
			);
create index on msmts (id);

-- Test SCADA device suite for ctown complete network
-- Control components
INSERT INTO channels (net_id_str, msmt_t) VALUES 
	('PU1', 'C'), 	('PU2', 'C'), 	('V2', 'C'), 	('PU4', 'C'), ('PU3', 'C'), 
	('PU5', 'C'), 	('PU6', 'C'), 	('PU7', 'C'), 	('PU8', 'C'),
	('PU10', 'C'),  	('PU11', 'C'),
    ('V2', 'V')
	;
-- Online sensors - tank levels and limits
INSERT INTO channels (net_id_str, msmt_t, l_lim, r_lim) VALUES
	('R1', 'L', 0, 0),  -- reservior
	('T1', 'L', 71.5, 78), ('T2', 'L', 65, 70.9), ('T3', 'L', 112.9, 112.9+6.75), 
	('T4', 'L', 132.5, 132.5+4.7), ('T5', 'L', 105.8, 105.8+4.5), 
	('T6', 'L', 101.5, 101.5+5.5), ('T7', 'L', 102, 102+5) ;


-- Online sensors - flowrates
INSERT INTO channels (net_id_str, msmt_t, l_lim, r_lim) VALUES
	('PU1', 'F', -200, 200), ('PU2', 'F', -200, 200), ('PU3', 'F', -200, 200), 
	('PU4', 'F', -200, 200), ('PU5', 'F', -200, 200), 
	('PU6', 'F', -200, 200), ('PU7', 'F', -200, 200), ('PU8', 'F', -200, 200), 
	('PU9', 'F', -200, 200), ('PU10', 'F', -200, 200), ('PU11', 'F', -200, 200),

	('P98', 'Q', -200, 200), ('P18', 'Q', -200, 200),('P946', 'Q', -200, 200),
	('P1016', 'Q', -200, 200), ('P138', 'Q', -200, 200),
	('P383', 'Q', -200, 200);

-- Online sensors - pressures
INSERT INTO channels (net_id_str, msmt_t) VALUES
	('J144', 'P'), ('J253', 'P'), ('J330', 'P'), ('J168', 'P'), ('J221', 'P'),
	('J28', 'P'),

    ('J273', 'B'),
    ('J269', 'B'),
    ('J274', 'B'),
    ('J292', 'B'),
    ('J256', 'B'),
    ('J415', 'B'),
    ('J291', 'B'),
    ('J304', 'B'),
    ('J306', 'B'),
    ('J317', 'B'),
    ('J323', 'B')
    ; 


	
