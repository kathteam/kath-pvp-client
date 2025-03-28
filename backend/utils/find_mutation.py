import mysql.connector

def get_mutations(chromosome, position):
    try:
        conn = mysql.connector.connect(
            host='172.20.128.1',
            user='root',
            password='password',
            database='test'
        )
        cursor = conn.cursor(dictionary=True)

        # Query to fetch mutations
        query = """
        (
            SELECT * FROM mutation_sorted
            WHERE Chromosome = %s
            AND Position = %s
            ORDER BY Position ASC
        );
        """
        
        cursor.execute(query, (chromosome, position))
        results = cursor.fetchall()

        if results:
            for row in results:
                print(row)
        else:
            print("No matching mutations found.")
        
    except mysql.connector.Error as err:
        print(f"Error: {err}")
    
    finally:
        if 'conn' in locals() and conn.is_connected():
            cursor.close()
            conn.close()



get_mutations(1, 4748753)