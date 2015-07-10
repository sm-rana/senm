set search_path to ctown_scada,public;
DROP TABLE IF EXISTS channels; 
-- DROP TABLE IF EXISTS Distros;

CREATE TABLE channels (  
			id		SERIAL, 
			net_id_str	VARCHAR(255) NOT NULL,
			msmt_t	CHAR NOT NULL,
			l_lim	DOUBLE PRECISION,
			r_lim	DOUBLE PRECISION,
			info	VARCHAR(255),
            stde    DOUBLE PRECISION,
			PRIMARY KEY (id)
			 
			); 
create index on channels (id);

-- Test SCADA device suite for ctown complete network

-- Control components
INSERT INTO channels (net_id_str, msmt_t) VALUES 
	('PU1', 'C'), 	('PU2', 'C'), 	('V2', 'C'), 	('PU4', 'C'), ('PU3', 'C'), 
	('PU5', 'C'), 	('PU6', 'C'), 	('PU7', 'C'), 	('PU8', 'C'),
	('PU10', 'C'),  	('PU11', 'C'),
    ('V2', 'V')
	;
-- Online sensors - tank levels and limits
INSERT INTO channels (net_id_str, msmt_t, l_lim, r_lim, stde) VALUES
	('R1', 'L', 0, 0, 1),  -- reservior
	('T1', 'L', 71.5, 78, 1), ('T2', 'L', 65, 70.9, 1), ('T3', 'L', 112.9,
        112.9+6.75, 1), 
	('T4', 'L', 132.5, 132.5+4.7, 1), ('T5', 'L', 105.8, 105.8+4.5, 1), 
	('T6', 'L', 101.5, 101.5+5.5, 1), ('T7', 'L', 102, 102+5, 1) ;


-- Online sensors - flowrates
INSERT INTO channels (net_id_str, msmt_t, l_lim, r_lim, stde) VALUES
	('PU1', 'F', -200, 200, 3), ('PU2', 'F', -200, 200, 3), ('PU3', 'F', -200,
        200, 3), 
	('PU4', 'F', -200, 200, 3), ('PU5', 'F', -200, 200, 3), 
	('PU6', 'F', -200, 200, 3), ('PU7', 'F', -200, 200, 3), ('PU8', 'F', -200, 200, 3), 
	('PU9', 'F', -200, 200, 3), ('PU10', 'F', -200, 200, 3), ('PU11', 'F', -200, 200, 3),

	('P98', 'Q', -200, 200, 3), ('P18', 'Q', -200, 200, 3),('P946', 'Q', -200, 200, 3),
	('P1016', 'Q', -200, 200, 3), ('P138', 'Q', -200, 200, 3),
	('P383', 'Q', -200, 200, 3);

-- Online sensors - pressures
INSERT INTO channels (net_id_str, msmt_t, l_lim, r_lim, stde) VALUES
	('J144', 'P', -200, 200, 5), ('J253', 'P', -200, 200, 5), ('J330', 'P', -200, 200, 5), ('J168', 'P', -200, 200, 5), ('J221', 'P', -200, 200, 5),
	('J28', 'P', -200, 200, 5),

    ('J273', 'B', -200, 200, 5),
    ('J269', 'B', -200, 200, 5),
    ('J274', 'B', -200, 200, 5),
    ('J292', 'B', -200, 200, 5),
    ('J256', 'B', -200, 200, 5),
    ('J415', 'B', -200, 200, 5),
    ('J291', 'B', -200, 200, 5),
    ('J304', 'B', -200, 200, 5),
    ('J306', 'B', -200, 200, 5),
    ('J317', 'B', -200, 200, 5),
    ('J323', 'B', -200, 200, 5)
    ; 


	
DROP TABLE IF EXISTS msmts; 
CREATE TABLE msmts (
			id		SERIAL, 
			time	TIMESTAMP NOT NULL, 
			value	DOUBLE PRECISION NOT NULL, 
			cid		INTEGER NOT NULL,
			PRIMARY KEY (id),
			CONSTRAINT fk_chnl FOREIGN KEY (cid) REFERENCES Channels(id) 
			);
create index on msmts (id);

