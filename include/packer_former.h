/**
 * @file packet_former.h
 * @author Riley Kitchenka (https://github.com/DrGrandmaster)
 * @brief Defines a class for forming SCHP5 packets
 * @version 0.1
 * @date 2025-10-02
 */
# ifndef PACKET_FORMER_H
# define PACKET_FORMER_H

#include "packet.h"
#include "ring.hpp"
#include "crcinterface.h"
#include <cstdint>
#include <queue>

/**
 * @brief Forms packets according to the SCHP5 protocol for packet radios
 */
class PacketFormer {

    private:

    /**
    * @brief The address of the sender
    */
    const uint32_t senderAddress;

    /**
    * @brief The address of the reciever
    */
    const uint32_t recieverAddress;

    /**
    * @brief The current length of the data portion of the packets
    */
    uint8_t dataLength;

    /**
     * @brief A pointer to the crc generator
     */
    crcutil_interface::CRC * crcGenny;

    public:

    /**
     * @brief Construct a new packet former
     */
    PacketFormer(
        uint32_t sender_address,
        uint32_t reciever_address,
        uint8_t data_length);


    /**
     * @brief Set Data Length
     */
    void setDataLength(uint8_t data_length);

    /**
     * @brief Form the next packet
     * @param data A pointer to the input data buffer
     * @param sequence_number The sequence number of the current packet
     * @return A well formed SCHP5 data packet
     */
    Packet formDataPacket(std::vector<uint8_t> * data, uint8_t sequence_number);

    /**
     * @brief Form a retransmitted packet
     * @param data A pointer to the input data buffer
     * @param sequence_number The sequence number of the current packet
     * @return A well formed SCHP5 retransmitted packet
     * 
     */
    Packet formRetransmitPacket(std::vector<uint8_t> * data, uint8_t sequence_num);

    /**
     * @brief Form a busy start packet
     * @return A well formed SCHP5 busy start packet
     */
    Packet formBusyStartPacket();

    /**
     * @brief Form a busy end packet
     * @return A well formed SCHP5 busy end packet
     */
    Packet formBusyEndPacket();

    /**
     * @brief Form a center frequency change packet
     * @param center_freq The center frequency to change to, as an IEEE 754 32-bit floating-point
     * @return A well formed SCHP5 center frequency change packet
     */
    Packet formCenterFrequencyPacket(const float center_freq);

    /**
     * @brief Form a modulation change packet
     * @param modulation_type The type of modulation to switch to, as specified in the SCHP5 protocol, usally log_2(QAM order - 1)
     * @return A well formed SCHP5 modulation change packet
     */
    Packet formModulationChangePacket(const uint8_t modulation_type);

    /**
     * @brief Form a transaction start packet
     * @param center_freq The center frequency to use, as an IEEE 754 32-bit floating-point
     * @param modulation_type The type of modulation to use, as specified in the SCHP5 protocol, usally log_2(QAM order - 1)
     * @return A well formed SCHP5 transaction start packet
     */
    Packet formTransactionStartPacket(const float center_freq, const uint8_t modulation_type);

    /**
     * @brief Form a transaction restart packet
     * @param center_freq The center frequency to use, as an IEEE 754 32-bit floating-point
     * @param modulation_type The type of modulation to use, as specified in the SCHP5 protocol, usally log_2(QAM order - 1)
     * @return A well formed SCHP5 transaction restart packet
     */
    Packet formTransactionRestartPacket(const float center_freq, const uint8_t modulation_type);

    /**
     * @brief Form a transaction transfer packet
     * @param center_freq The center frequency to use, as an IEEE 754 32-bit floating-point
     * @param modulation_type The type of modulation to use, as specified in the SCHP5 protocol, usally log_2(QAM order - 1)
     * @return A well formed SCHP5 transaction transfer packet
     */
    Packet formTransactionTransferPacket(const float center_freq, const uint8_t modulation_type);

    /**
     * @brief Form a transaction dropped packet
     * @return A well formed SCHP5 transaction dropped packet
     */
    Packet formTransactionDroppedPacket();

    /**
     * @brief Form a transaction end packet
     * @return A well formed SCHP5 transaction end packet
     */
    Packet formTransactionEndPacket();

};

#endif